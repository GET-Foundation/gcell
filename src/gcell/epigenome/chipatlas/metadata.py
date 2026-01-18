"""ChIP-Atlas metadata management with SQLite caching.

This module handles downloading, caching, and querying ChIP-Atlas metadata
using SQLite for fast searches.
"""

import io
import sqlite3
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import TYPE_CHECKING

import pandas as pd

from gcell import _settings

if TYPE_CHECKING:
    from collections.abc import Iterator

# ChIP-Atlas metadata URLs
BASE_URL = "https://chip-atlas.dbcls.jp/data"
METADATA_BASE_URL = f"{BASE_URL}/metadata"

METADATA_FILES = {
    "experimentList": f"{METADATA_BASE_URL}/experimentList.tab",
    "fileList": f"{METADATA_BASE_URL}/fileList.tab",
    "antigenList": f"{METADATA_BASE_URL}/antigenList.tab",
    "celltypeList": f"{METADATA_BASE_URL}/celltypeList.tab",
    "analysisList": f"{METADATA_BASE_URL}/analysisList.tab",
}

# Column names for each metadata file
EXPERIMENT_COLUMNS = [
    "experiment_id",
    "assembly",
    "antigen_class",
    "antigen",
    "cell_type_class",
    "cell_type",
    "cell_type_description",
    "processing_logs",
    "title",
    "metadata",
]

FILE_COLUMNS = [
    "file_name",
    "assembly",
    "antigen_class",
    "antigen",
    "cell_type_class",
    "cell_type",
    "threshold",
    "experiment_ids",
]

ANTIGEN_COLUMNS = [
    "assembly",
    "antigen_class",
    "antigen",
    "experiment_count",
    "experiment_ids",
]

CELLTYPE_COLUMNS = [
    "assembly",
    "cell_type_class",
    "cell_type",
    "experiment_count",
    "experiment_ids",
]

ANALYSIS_COLUMNS = [
    "antigen",
    "cell_type_classes",
    "has_target_genes",
    "assembly",
]

# Supported genome assemblies
SUPPORTED_ASSEMBLIES = [
    "hg38",
    "hg19",
    "mm10",
    "mm9",
    "rn6",
    "dm6",
    "dm3",
    "ce11",
    "ce10",
    "sacCer3",
]

# Antigen classes
ANTIGEN_CLASSES = [
    "TFs and others",
    "Histone",
    "Input control",
    "RNA polymerase",
    "ATAC-Seq",
    "DNase-seq",
    "Bisulfite-Seq",
]


def get_chipatlas_dir() -> Path:
    """Get the ChIP-Atlas data directory."""
    chipatlas_dir = Path(_settings.get_setting("annotation_dir")) / "chipatlas"
    chipatlas_dir.mkdir(parents=True, exist_ok=True)
    return chipatlas_dir


def get_db_path() -> Path:
    """Get path to the SQLite database."""
    return get_chipatlas_dir() / "chipatlas.db"


# Default row limits for lite mode (enough for exploration)
LITE_MODE_ROWS = 5000


class ChipAtlasMetadata:
    """Manages ChIP-Atlas metadata with SQLite caching.

    This class downloads and caches ChIP-Atlas metadata in a SQLite database
    for efficient querying. The database is created on first access and
    can be refreshed manually.

    Parameters
    ----------
    force_refresh : bool, optional
        If True, force re-download of metadata files. Default is False.
    max_workers : int, optional
        Maximum number of parallel download workers. Default is 4.
    metadata_mode : {"full", "lite"}, optional
        Mode for metadata caching:
        - "full": Download complete metadata (~500MB). Slower first-time setup
          but comprehensive coverage of all experiments.
        - "lite": Download partial metadata (~1MB). Quick setup for exploration
          but limited to ~5000 experiments per table.
        Default is "full".

    Attributes
    ----------
    db_path : Path
        Path to the SQLite database file.

    Examples
    --------
    >>> # Full mode (default) - complete metadata
    >>> meta = ChipAtlasMetadata()

    >>> # Lite mode - quick setup for exploration
    >>> meta = ChipAtlasMetadata(metadata_mode="lite")

    >>> # Get all available antigens for human
    >>> antigens = meta.get_antigens(assembly="hg38")

    >>> # Upgrade from lite to full mode
    >>> meta.upgrade_to_full()
    """

    def __init__(
        self,
        force_refresh: bool = False,
        max_workers: int = 4,
        metadata_mode: str = "full",
        max_rows: int | None = None,
    ):
        if metadata_mode not in ("full", "lite"):
            raise ValueError("metadata_mode must be 'full' or 'lite'")

        self.db_path = get_db_path()
        self.max_workers = max_workers
        self.metadata_mode = metadata_mode
        # max_rows overrides metadata_mode for testing
        self.max_rows = max_rows

        # Check if we need to initialize or upgrade
        if force_refresh:
            self._initialize_db()
        elif not self._is_initialized():
            if self._should_upgrade():
                # User requested full mode but database is lite
                print("Upgrading metadata from lite to full mode...")
            self._initialize_db()

    def _db_exists(self) -> bool:
        """Check if database exists and has all required tables."""
        if not self.db_path.exists():
            return False

        try:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.execute(
                    "SELECT name FROM sqlite_master WHERE type='table'"
                )
                tables = {row[0] for row in cursor.fetchall()}
                required = {
                    "experiments",
                    "files",
                    "antigens",
                    "celltypes",
                    "metadata_info",
                }
                return required.issubset(tables)
        except sqlite3.Error:
            return False

    def _should_upgrade(self) -> bool:
        """Check if we should upgrade from lite to full mode.

        Returns True only if user requested full mode AND current DB is lite.
        """
        if self.metadata_mode != "full":
            return False

        # Check current mode in database
        current_mode = self._get_stored_mode()
        # Only upgrade if explicitly lite mode (not if already full or unknown)
        return current_mode == "lite"

    def _is_initialized(self) -> bool:
        """Check if the database is properly initialized for the requested mode.

        Returns True if:
        - DB exists with all required tables
        - Mode matches requested mode (or requested is lite and stored is full)
        """
        if not self._db_exists():
            return False

        stored_mode = self._get_stored_mode()
        if stored_mode is None:
            return False

        # If requesting lite, full mode is also acceptable
        if self.metadata_mode == "lite":
            return stored_mode in ("lite", "full")

        # If requesting full, only full is acceptable
        return stored_mode == self.metadata_mode

    def _get_stored_mode(self) -> str | None:
        """Get the metadata mode stored in the database."""
        if not self.db_path.exists():
            return None

        try:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.execute(
                    "SELECT value FROM metadata_info WHERE key = 'metadata_mode'"
                )
                row = cursor.fetchone()
                return row[0] if row else None
        except sqlite3.Error:
            return None

    def _initialize_db(self) -> None:
        """Initialize the SQLite database with metadata."""
        effective_mode = self.metadata_mode
        if self.max_rows is not None:
            effective_mode = "custom"

        if effective_mode == "lite":
            print("Initializing ChIP-Atlas metadata database (lite mode)...")
            print("Downloading ~1MB of metadata for quick exploration.")
            print("Use upgrade_to_full() or metadata_mode='full' for complete data.")
        elif effective_mode == "full":
            print("Initializing ChIP-Atlas metadata database (full mode)...")
            print("Downloading complete metadata (~500MB). This may take a few minutes.")
        else:
            print("Initializing ChIP-Atlas metadata database...")

        # Download all metadata files in parallel
        metadata = self._download_all_metadata()

        # Create database
        with sqlite3.connect(self.db_path) as conn:
            self._create_tables(conn)
            self._populate_tables(conn, metadata)
            self._create_indexes(conn)
            # Store metadata mode
            conn.execute(
                "INSERT OR REPLACE INTO metadata_info VALUES (?, ?)",
                ("metadata_mode", effective_mode),
            )

        print("ChIP-Atlas metadata database initialized successfully.")

    def _download_metadata_file(
        self,
        name: str,
        url: str,
        max_bytes: int | None = None,
        max_rows: int | None = None,
    ) -> tuple[str, pd.DataFrame]:
        """Download and parse a single metadata file.

        Parameters
        ----------
        name : str
            Name of the metadata file
        url : str
            URL to download from
        max_bytes : int, optional
            Maximum bytes to download (for streaming partial data).
            If None, downloads the full file.
        max_rows : int, optional
            Maximum rows to keep after parsing.
        """
        import urllib.request

        print(f"  Downloading {name}...")

        try:
            req = urllib.request.Request(url)
            if max_bytes is not None:
                # Use HTTP Range header for partial download
                req.add_header("Range", f"bytes=0-{max_bytes}")

            with urllib.request.urlopen(req, timeout=60) as response:
                content = response.read()
        except Exception as e:
            raise RuntimeError(f"Failed to download {name}: {e}") from e

        # If partial download, remove potentially incomplete last line
        if max_bytes is not None:
            content_str = content.decode("utf-8", errors="ignore")
            lines = content_str.split("\n")
            # Remove last line if it might be incomplete
            if lines and not content_str.endswith("\n"):
                lines = lines[:-1]
            content = "\n".join(lines).encode("utf-8")

        # Column definitions for each file type
        columns_map = {
            "experimentList": EXPERIMENT_COLUMNS,
            "fileList": FILE_COLUMNS,
            "antigenList": ANTIGEN_COLUMNS,
            "celltypeList": CELLTYPE_COLUMNS,
            "analysisList": ANALYSIS_COLUMNS,
        }

        if name not in columns_map:
            raise ValueError(f"Unknown metadata file: {name}")

        # Parse with on_bad_lines='skip' to handle potential issues
        # Use index_col=False to prevent pandas from using extra columns as index
        # Use usecols to only read the columns we need (files may have extra fields)
        num_cols = len(columns_map[name])
        df = pd.read_csv(
            io.BytesIO(content),
            sep="\t",
            names=columns_map[name],
            usecols=range(num_cols),
            index_col=False,
            dtype=str,
            na_values=["", "-"],
            low_memory=False,
            on_bad_lines="skip",
        )

        # Apply max_rows limit if specified
        if max_rows is not None and len(df) > max_rows:
            df = df.head(max_rows)

        return name, df

    def _download_all_metadata(self) -> dict[str, pd.DataFrame]:
        """Download all metadata files in parallel."""
        metadata = {}

        # Determine row limit based on mode
        max_rows = self.max_rows  # Explicit max_rows takes precedence
        if max_rows is None and self.metadata_mode == "lite":
            max_rows = LITE_MODE_ROWS

        # Calculate max_bytes if max_rows is set
        # Estimate ~500 bytes per row for experiment data
        max_bytes = None
        if max_rows is not None:
            max_bytes = max_rows * 500

        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(
                    self._download_metadata_file, name, url, max_bytes, max_rows
                ): name
                for name, url in METADATA_FILES.items()
            }

            for future in as_completed(futures):
                name, df = future.result()
                metadata[name] = df

        return metadata

    def _create_tables(self, conn: sqlite3.Connection) -> None:
        """Create database tables."""
        conn.executescript(
            """
            DROP TABLE IF EXISTS experiments;
            DROP TABLE IF EXISTS files;
            DROP TABLE IF EXISTS antigens;
            DROP TABLE IF EXISTS celltypes;
            DROP TABLE IF EXISTS metadata_info;

            CREATE TABLE experiments (
                experiment_id TEXT PRIMARY KEY,
                assembly TEXT NOT NULL,
                antigen_class TEXT,
                antigen TEXT,
                cell_type_class TEXT,
                cell_type TEXT,
                cell_type_description TEXT,
                processing_logs TEXT,
                title TEXT,
                metadata TEXT
            );

            CREATE TABLE files (
                file_name TEXT PRIMARY KEY,
                assembly TEXT NOT NULL,
                antigen_class TEXT,
                antigen TEXT,
                cell_type_class TEXT,
                cell_type TEXT,
                threshold INTEGER,
                experiment_ids TEXT
            );

            CREATE TABLE antigens (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                assembly TEXT NOT NULL,
                antigen_class TEXT,
                antigen TEXT,
                experiment_count INTEGER,
                experiment_ids TEXT
            );

            CREATE TABLE celltypes (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                assembly TEXT NOT NULL,
                cell_type_class TEXT,
                cell_type TEXT,
                experiment_count INTEGER,
                experiment_ids TEXT
            );

            CREATE TABLE metadata_info (
                key TEXT PRIMARY KEY,
                value TEXT
            );
        """
        )

    def _populate_tables(
        self, conn: sqlite3.Connection, metadata: dict[str, pd.DataFrame]
    ) -> None:
        """Populate database tables with metadata."""
        # Experiments
        if "experimentList" in metadata:
            df = metadata["experimentList"].dropna(subset=["experiment_id"])
            df.to_sql("experiments", conn, if_exists="replace", index=False)

        # Files
        if "fileList" in metadata:
            df = metadata["fileList"].dropna(subset=["file_name"])
            # Convert threshold to integer
            df["threshold"] = pd.to_numeric(df["threshold"], errors="coerce")
            df.to_sql("files", conn, if_exists="replace", index=False)

        # Antigens
        if "antigenList" in metadata:
            df = metadata["antigenList"].dropna(subset=["antigen"]).copy()
            df["experiment_count"] = pd.to_numeric(
                df["experiment_count"], errors="coerce"
            )
            df.to_sql("antigens", conn, if_exists="replace", index=False)

        # Cell types
        if "celltypeList" in metadata:
            df = metadata["celltypeList"].dropna(subset=["cell_type"]).copy()
            df["experiment_count"] = pd.to_numeric(
                df["experiment_count"], errors="coerce"
            )
            df.to_sql("celltypes", conn, if_exists="replace", index=False)

        # Store metadata info
        import datetime

        conn.execute(
            "INSERT OR REPLACE INTO metadata_info VALUES (?, ?)",
            ("last_updated", datetime.datetime.now(datetime.UTC).isoformat()),
        )

    def _create_indexes(self, conn: sqlite3.Connection) -> None:
        """Create indexes for fast querying."""
        conn.executescript(
            """
            CREATE INDEX IF NOT EXISTS idx_exp_assembly ON experiments(assembly);
            CREATE INDEX IF NOT EXISTS idx_exp_antigen ON experiments(antigen);
            CREATE INDEX IF NOT EXISTS idx_exp_antigen_class ON experiments(antigen_class);
            CREATE INDEX IF NOT EXISTS idx_exp_celltype ON experiments(cell_type);
            CREATE INDEX IF NOT EXISTS idx_exp_celltype_class ON experiments(cell_type_class);
            CREATE INDEX IF NOT EXISTS idx_exp_composite
                ON experiments(assembly, antigen_class, cell_type_class);

            CREATE INDEX IF NOT EXISTS idx_files_assembly ON files(assembly);
            CREATE INDEX IF NOT EXISTS idx_files_antigen ON files(antigen);
            CREATE INDEX IF NOT EXISTS idx_files_celltype ON files(cell_type);
            CREATE INDEX IF NOT EXISTS idx_files_threshold ON files(threshold);

            CREATE INDEX IF NOT EXISTS idx_antigens_assembly ON antigens(assembly);
            CREATE INDEX IF NOT EXISTS idx_antigens_class ON antigens(antigen_class);

            CREATE INDEX IF NOT EXISTS idx_celltypes_assembly ON celltypes(assembly);
            CREATE INDEX IF NOT EXISTS idx_celltypes_class ON celltypes(cell_type_class);
        """
        )

    def _execute_query(self, query: str, params: tuple = ()) -> "Iterator[sqlite3.Row]":
        """Execute a query and return results."""
        with sqlite3.connect(self.db_path) as conn:
            conn.row_factory = sqlite3.Row
            cursor = conn.execute(query, params)
            return cursor.fetchall()

    def search_experiments(
        self,
        assembly: str | None = None,
        antigen: str | None = None,
        antigen_class: str | None = None,
        cell_type: str | None = None,
        cell_type_class: str | None = None,
        title_contains: str | None = None,
        limit: int | None = None,
        sort_by_library_size: bool = False,
    ) -> pd.DataFrame:
        """Search for experiments matching criteria.

        Parameters
        ----------
        assembly : str, optional
            Genome assembly (hg38, hg19, mm10, etc.)
        antigen : str, optional
            Antigen/target name (e.g., TP53, H3K27ac)
        antigen_class : str, optional
            Antigen class (TFs and others, Histone, etc.)
        cell_type : str, optional
            Cell type name
        cell_type_class : str, optional
            Cell type class (Blood, Breast, etc.)
        title_contains : str, optional
            Search for experiments with title containing this string
        limit : int, optional
            Maximum number of results to return
        sort_by_library_size : bool, optional
            If True, sort results by library size (read count) descending.
            Library size is extracted from processing_logs field.

        Returns
        -------
        pd.DataFrame
            DataFrame with matching experiments
        """
        conditions = []
        params = []

        if assembly:
            conditions.append("assembly = ?")
            params.append(assembly)
        if antigen:
            conditions.append("antigen LIKE ?")
            params.append(f"%{antigen}%")
        if antigen_class:
            conditions.append("antigen_class = ?")
            params.append(antigen_class)
        if cell_type:
            conditions.append("cell_type LIKE ?")
            params.append(f"%{cell_type}%")
        if cell_type_class:
            conditions.append("cell_type_class = ?")
            params.append(cell_type_class)
        if title_contains:
            conditions.append("title LIKE ?")
            params.append(f"%{title_contains}%")

        query = "SELECT * FROM experiments"
        if conditions:
            query += " WHERE " + " AND ".join(conditions)

        rows = self._execute_query(query, tuple(params))
        df = pd.DataFrame([dict(row) for row in rows])

        # Sort by library size if requested
        if sort_by_library_size and len(df) > 0 and "processing_logs" in df.columns:
            df["_library_size"] = df["processing_logs"].apply(_extract_library_size)
            df = df.sort_values("_library_size", ascending=False)
            df = df.drop(columns=["_library_size"])

        if limit:
            df = df.head(limit)

        return df


    def get_antigens(
        self,
        assembly: str | None = None,
        antigen_class: str | None = None,
        min_experiments: int = 0,
    ) -> pd.DataFrame:
        """Get available antigens/targets.

        Parameters
        ----------
        assembly : str, optional
            Filter by genome assembly
        antigen_class : str, optional
            Filter by antigen class
        min_experiments : int, optional
            Minimum number of experiments required

        Returns
        -------
        pd.DataFrame
            DataFrame with antigen information
        """
        conditions = ["experiment_count >= ?"]
        params = [min_experiments]

        if assembly:
            conditions.append("assembly = ?")
            params.append(assembly)
        if antigen_class:
            conditions.append("antigen_class = ?")
            params.append(antigen_class)

        query = f"""
            SELECT DISTINCT assembly, antigen_class, antigen, experiment_count
            FROM antigens
            WHERE {" AND ".join(conditions)}
            ORDER BY experiment_count DESC
        """

        rows = self._execute_query(query, tuple(params))
        return pd.DataFrame([dict(row) for row in rows])

    def get_celltypes(
        self,
        assembly: str | None = None,
        cell_type_class: str | None = None,
        min_experiments: int = 0,
    ) -> pd.DataFrame:
        """Get available cell types.

        Parameters
        ----------
        assembly : str, optional
            Filter by genome assembly
        cell_type_class : str, optional
            Filter by cell type class
        min_experiments : int, optional
            Minimum number of experiments required

        Returns
        -------
        pd.DataFrame
            DataFrame with cell type information
        """
        conditions = ["experiment_count >= ?"]
        params = [min_experiments]

        if assembly:
            conditions.append("assembly = ?")
            params.append(assembly)
        if cell_type_class:
            conditions.append("cell_type_class = ?")
            params.append(cell_type_class)

        query = f"""
            SELECT DISTINCT assembly, cell_type_class, cell_type, experiment_count
            FROM celltypes
            WHERE {" AND ".join(conditions)}
            ORDER BY experiment_count DESC
        """

        rows = self._execute_query(query, tuple(params))
        return pd.DataFrame([dict(row) for row in rows])

    def get_celltype_classes(self, assembly: str | None = None) -> list[str]:
        """Get available cell type classes.

        Parameters
        ----------
        assembly : str, optional
            Filter by genome assembly

        Returns
        -------
        list[str]
            List of cell type class names
        """
        query = "SELECT DISTINCT cell_type_class FROM celltypes"
        params = []

        if assembly:
            query += " WHERE assembly = ?"
            params.append(assembly)

        query += " ORDER BY cell_type_class"

        rows = self._execute_query(query, tuple(params))
        return [row["cell_type_class"] for row in rows if row["cell_type_class"]]

    def get_antigen_classes(self, assembly: str | None = None) -> list[str]:
        """Get available antigen classes.

        Parameters
        ----------
        assembly : str, optional
            Filter by genome assembly

        Returns
        -------
        list[str]
            List of antigen class names
        """
        query = "SELECT DISTINCT antigen_class FROM antigens"
        params = []

        if assembly:
            query += " WHERE assembly = ?"
            params.append(assembly)

        query += " ORDER BY antigen_class"

        rows = self._execute_query(query, tuple(params))
        return [row["antigen_class"] for row in rows if row["antigen_class"]]

    def get_assemblies(self) -> list[str]:
        """Get available genome assemblies.

        Returns
        -------
        list[str]
            List of genome assembly names
        """
        rows = self._execute_query(
            "SELECT DISTINCT assembly FROM experiments ORDER BY assembly"
        )
        return [row["assembly"] for row in rows if row["assembly"]]

    def get_experiment_count(
        self,
        assembly: str | None = None,
        antigen_class: str | None = None,
        cell_type_class: str | None = None,
    ) -> int:
        """Get count of experiments matching criteria.

        Parameters
        ----------
        assembly : str, optional
            Filter by genome assembly
        antigen_class : str, optional
            Filter by antigen class
        cell_type_class : str, optional
            Filter by cell type class

        Returns
        -------
        int
            Number of matching experiments
        """
        conditions = []
        params = []

        if assembly:
            conditions.append("assembly = ?")
            params.append(assembly)
        if antigen_class:
            conditions.append("antigen_class = ?")
            params.append(antigen_class)
        if cell_type_class:
            conditions.append("cell_type_class = ?")
            params.append(cell_type_class)

        query = "SELECT COUNT(*) as count FROM experiments"
        if conditions:
            query += " WHERE " + " AND ".join(conditions)

        rows = self._execute_query(query, tuple(params))
        return rows[0]["count"]

    def get_last_updated(self) -> str | None:
        """Get the last update time of the metadata database.

        Returns
        -------
        str | None
            ISO format datetime string or None if not available
        """
        rows = self._execute_query(
            "SELECT value FROM metadata_info WHERE key = 'last_updated'"
        )
        return rows[0]["value"] if rows else None

    def refresh(self) -> None:
        """Force refresh of the metadata database."""
        self._initialize_db()

    def upgrade_to_full(self) -> None:
        """Upgrade from lite mode to full metadata.

        Downloads the complete metadata if currently in lite mode.
        """
        current_mode = self._get_stored_mode()
        if current_mode == "full":
            print("Metadata is already in full mode.")
            return

        print("Upgrading metadata from lite to full mode...")
        self.metadata_mode = "full"
        self._initialize_db()

    def get_mode(self) -> str:
        """Get the current metadata mode.

        Returns
        -------
        str
            Current mode: "full", "lite", or "unknown"
        """
        mode = self._get_stored_mode()
        return mode if mode else "unknown"

    def is_lite_mode(self) -> bool:
        """Check if database is in lite mode.

        Returns
        -------
        bool
            True if in lite mode, False otherwise
        """
        return self._get_stored_mode() == "lite"


def _extract_library_size(processing_logs: str | None) -> int:
    """Extract library size (read count) from processing_logs field.

    The processing_logs field format is: read_count,alignment_rate,dup_rate,...

    Parameters
    ----------
    processing_logs : str or None
        The processing_logs field value

    Returns
    -------
    int
        Library size (read count), or 0 if not parseable
    """
    if not processing_logs or not isinstance(processing_logs, str):
        return 0
    try:
        # First value is the read count
        first_val = processing_logs.split(",")[0]
        return int(first_val)
    except (ValueError, IndexError):
        return 0
