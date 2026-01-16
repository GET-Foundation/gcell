"""Tests for the ChIP-Atlas query module."""

import sqlite3
from unittest.mock import patch

import pandas as pd
import pytest

from gcell.epigenome.chipatlas.metadata import (
    ANTIGEN_CLASSES,
    SUPPORTED_ASSEMBLIES,
    ChipAtlasMetadata,
)
from gcell.epigenome.chipatlas.query import (
    THRESHOLD_MAP,
    ChipAtlas,
    ChipAtlasExperiment,
)


class TestChipAtlasExperiment:
    """Tests for ChipAtlasExperiment dataclass."""

    def test_experiment_creation(self):
        """Test creating an experiment object."""
        exp = ChipAtlasExperiment(
            experiment_id="SRX123456",
            assembly="hg38",
            antigen="TP53",
            antigen_class="TFs and others",
            cell_type="HeLa",
            cell_type_class="Breast",
            title="TP53 ChIP-seq in HeLa cells",
        )

        assert exp.experiment_id == "SRX123456"
        assert exp.assembly == "hg38"
        assert exp.antigen == "TP53"

    def test_bigwig_url(self):
        """Test BigWig URL generation."""
        exp = ChipAtlasExperiment(
            experiment_id="SRX123456",
            assembly="hg38",
        )

        url = exp.bigwig_url
        assert "SRX123456.bw" in url
        assert "hg38" in url

    def test_peak_url(self):
        """Test peak URL generation."""
        exp = ChipAtlasExperiment(
            experiment_id="SRX123456",
            assembly="hg38",
        )

        # Test different thresholds
        for threshold, thresh_str in THRESHOLD_MAP.items():
            url = exp.peak_url(threshold=threshold)
            assert f"{thresh_str}.bed" in url
            assert "SRX123456" in url

    def test_bigbed_url(self):
        """Test BigBed URL generation."""
        exp = ChipAtlasExperiment(
            experiment_id="SRX123456",
            assembly="hg38",
        )

        url = exp.bigbed_url(threshold=10)
        assert "SRX123456.10.bb" in url


class TestChipAtlasMetadataConstants:
    """Tests for metadata constants."""

    def test_supported_assemblies(self):
        """Test that common assemblies are supported."""
        assert "hg38" in SUPPORTED_ASSEMBLIES
        assert "hg19" in SUPPORTED_ASSEMBLIES
        assert "mm10" in SUPPORTED_ASSEMBLIES
        assert "mm9" in SUPPORTED_ASSEMBLIES

    def test_antigen_classes(self):
        """Test that expected antigen classes are defined."""
        assert "TFs and others" in ANTIGEN_CLASSES
        assert "Histone" in ANTIGEN_CLASSES
        assert "ATAC-Seq" in ANTIGEN_CLASSES
        assert "DNase-seq" in ANTIGEN_CLASSES


class TestChipAtlasMetadataOffline:
    """Offline tests for ChipAtlasMetadata using mock database."""

    @pytest.fixture
    def mock_db(self, tmp_path):
        """Create a mock SQLite database for testing."""
        db_path = tmp_path / "test_chipatlas.db"

        with sqlite3.connect(db_path) as conn:
            # Create tables
            conn.executescript(
                """
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

            # Insert test data
            conn.executemany(
                """
                INSERT INTO experiments VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        "SRX000001",
                        "hg38",
                        "TFs and others",
                        "TP53",
                        "Breast",
                        "HeLa",
                        "HeLa cells",
                        "OK",
                        "TP53 ChIP-seq",
                        "{}",
                    ),
                    (
                        "SRX000002",
                        "hg38",
                        "Histone",
                        "H3K27ac",
                        "Blood",
                        "K562",
                        "K562 cells",
                        "OK",
                        "H3K27ac ChIP-seq",
                        "{}",
                    ),
                    (
                        "SRX000003",
                        "mm10",
                        "TFs and others",
                        "CTCF",
                        "Embryo",
                        "mESC",
                        "Mouse ESC",
                        "OK",
                        "CTCF ChIP-seq",
                        "{}",
                    ),
                ],
            )

            conn.executemany(
                """
                INSERT INTO antigens (assembly, antigen_class, antigen, experiment_count, experiment_ids)
                VALUES (?, ?, ?, ?, ?)
                """,
                [
                    ("hg38", "TFs and others", "TP53", 100, "SRX000001"),
                    ("hg38", "Histone", "H3K27ac", 500, "SRX000002"),
                    ("mm10", "TFs and others", "CTCF", 200, "SRX000003"),
                ],
            )

            conn.executemany(
                """
                INSERT INTO celltypes (assembly, cell_type_class, cell_type, experiment_count, experiment_ids)
                VALUES (?, ?, ?, ?, ?)
                """,
                [
                    ("hg38", "Breast", "HeLa", 50, "SRX000001"),
                    ("hg38", "Blood", "K562", 100, "SRX000002"),
                    ("mm10", "Embryo", "mESC", 75, "SRX000003"),
                ],
            )

            conn.execute(
                "INSERT INTO metadata_info VALUES (?, ?)",
                ("last_updated", "2024-01-01T00:00:00Z"),
            )

            # Create indexes
            conn.executescript(
                """
                CREATE INDEX idx_exp_assembly ON experiments(assembly);
                CREATE INDEX idx_exp_antigen ON experiments(antigen);
                CREATE INDEX idx_antigens_assembly ON antigens(assembly);
                CREATE INDEX idx_celltypes_assembly ON celltypes(assembly);
            """
            )

        return db_path

    def test_search_experiments_by_antigen(self, mock_db):
        """Test searching experiments by antigen."""
        with patch(
            "gcell.epigenome.chipatlas.metadata.get_db_path", return_value=mock_db
        ):
            meta = ChipAtlasMetadata.__new__(ChipAtlasMetadata)
            meta.db_path = mock_db
            meta.max_workers = 4

            results = meta.search_experiments(antigen="TP53")
            assert len(results) == 1
            assert results.iloc[0]["experiment_id"] == "SRX000001"

    def test_search_experiments_by_assembly(self, mock_db):
        """Test searching experiments by assembly."""
        with patch(
            "gcell.epigenome.chipatlas.metadata.get_db_path", return_value=mock_db
        ):
            meta = ChipAtlasMetadata.__new__(ChipAtlasMetadata)
            meta.db_path = mock_db
            meta.max_workers = 4

            results = meta.search_experiments(assembly="hg38")
            assert len(results) == 2

            results = meta.search_experiments(assembly="mm10")
            assert len(results) == 1

    def test_search_experiments_combined(self, mock_db):
        """Test searching with multiple criteria."""
        with patch(
            "gcell.epigenome.chipatlas.metadata.get_db_path", return_value=mock_db
        ):
            meta = ChipAtlasMetadata.__new__(ChipAtlasMetadata)
            meta.db_path = mock_db
            meta.max_workers = 4

            results = meta.search_experiments(assembly="hg38", antigen_class="Histone")
            assert len(results) == 1
            assert results.iloc[0]["antigen"] == "H3K27ac"

    def test_get_antigens(self, mock_db):
        """Test getting available antigens."""
        with patch(
            "gcell.epigenome.chipatlas.metadata.get_db_path", return_value=mock_db
        ):
            meta = ChipAtlasMetadata.__new__(ChipAtlasMetadata)
            meta.db_path = mock_db
            meta.max_workers = 4

            antigens = meta.get_antigens(assembly="hg38")
            assert len(antigens) == 2
            assert "TP53" in antigens["antigen"].values
            assert "H3K27ac" in antigens["antigen"].values

    def test_get_celltypes(self, mock_db):
        """Test getting available cell types."""
        with patch(
            "gcell.epigenome.chipatlas.metadata.get_db_path", return_value=mock_db
        ):
            meta = ChipAtlasMetadata.__new__(ChipAtlasMetadata)
            meta.db_path = mock_db
            meta.max_workers = 4

            celltypes = meta.get_celltypes(assembly="hg38")
            assert len(celltypes) == 2
            assert "HeLa" in celltypes["cell_type"].values
            assert "K562" in celltypes["cell_type"].values

    def test_get_celltype_classes(self, mock_db):
        """Test getting cell type classes."""
        with patch(
            "gcell.epigenome.chipatlas.metadata.get_db_path", return_value=mock_db
        ):
            meta = ChipAtlasMetadata.__new__(ChipAtlasMetadata)
            meta.db_path = mock_db
            meta.max_workers = 4

            classes = meta.get_celltype_classes(assembly="hg38")
            assert "Blood" in classes
            assert "Breast" in classes

    def test_get_experiment_count(self, mock_db):
        """Test getting experiment counts."""
        with patch(
            "gcell.epigenome.chipatlas.metadata.get_db_path", return_value=mock_db
        ):
            meta = ChipAtlasMetadata.__new__(ChipAtlasMetadata)
            meta.db_path = mock_db
            meta.max_workers = 4

            count = meta.get_experiment_count(assembly="hg38")
            assert count == 2

            count = meta.get_experiment_count(assembly="mm10")
            assert count == 1


class TestChipAtlasOffline:
    """Offline tests for ChipAtlas main class."""

    @pytest.fixture
    def mock_chipatlas(self, tmp_path):
        """Create a ChipAtlas instance with mock database."""
        db_path = tmp_path / "test_chipatlas.db"

        # Create mock database
        with sqlite3.connect(db_path) as conn:
            conn.executescript(
                """
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

            conn.execute(
                """
                INSERT INTO experiments VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    "SRX000001",
                    "hg38",
                    "TFs and others",
                    "TP53",
                    "Breast",
                    "HeLa",
                    "HeLa cells",
                    "OK",
                    "TP53 ChIP-seq",
                    "{}",
                ),
            )

            conn.execute(
                "INSERT INTO metadata_info VALUES (?, ?)",
                ("last_updated", "2024-01-01T00:00:00Z"),
            )

        # Create ChipAtlas with mock - we only need to patch in metadata module
        # since query.py imports get_db_path from metadata at call time
        with (
            patch(
                "gcell.epigenome.chipatlas.metadata.get_db_path", return_value=db_path
            ),
            patch.object(ChipAtlasMetadata, "_db_exists", return_value=True),
        ):
            ca = ChipAtlas(cache_dir=tmp_path)
            ca.metadata.db_path = db_path
            return ca

    def test_supported_assemblies_property(self, mock_chipatlas):
        """Test supported_assemblies property."""
        assemblies = mock_chipatlas.supported_assemblies
        assert isinstance(assemblies, list)
        assert "hg38" in assemblies

    def test_antigen_classes_property(self, mock_chipatlas):
        """Test antigen_classes property."""
        classes = mock_chipatlas.antigen_classes
        assert isinstance(classes, list)
        assert "TFs and others" in classes
        assert "Histone" in classes

    def test_search_returns_dataframe(self, mock_chipatlas):
        """Test that search returns a DataFrame."""
        results = mock_chipatlas.search(antigen="TP53")
        assert isinstance(results, pd.DataFrame)

    def test_search_as_experiments(self, mock_chipatlas):
        """Test search with as_experiments=True."""
        results = mock_chipatlas.search(antigen="TP53", as_experiments=True)
        assert isinstance(results, list)
        if len(results) > 0:
            assert isinstance(results[0], ChipAtlasExperiment)


class TestThresholdMap:
    """Tests for threshold mapping."""

    def test_threshold_values(self):
        """Test that threshold map contains expected values."""
        assert THRESHOLD_MAP[5] == "05"
        assert THRESHOLD_MAP[10] == "10"
        assert THRESHOLD_MAP[20] == "20"

    def test_threshold_keys(self):
        """Test threshold keys."""
        assert set(THRESHOLD_MAP.keys()) == {5, 10, 20}


class TestMetadataMode:
    """Tests for metadata mode functionality."""

    def test_invalid_metadata_mode(self):
        """Test that invalid metadata_mode raises ValueError."""
        with pytest.raises(ValueError, match="metadata_mode must be"):
            ChipAtlasMetadata(metadata_mode="invalid")

    def test_metadata_mode_stored_in_db(self, tmp_path):
        """Test that metadata mode is stored in the database."""
        db_path = tmp_path / "test_chipatlas.db"

        # Create mock database with mode stored
        with sqlite3.connect(db_path) as conn:
            conn.executescript(
                """
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

            # Store mode as "lite"
            conn.execute(
                "INSERT INTO metadata_info VALUES (?, ?)",
                ("metadata_mode", "lite"),
            )
            conn.execute(
                "INSERT INTO metadata_info VALUES (?, ?)",
                ("last_updated", "2024-01-01T00:00:00Z"),
            )

        # Test reading mode from database
        with patch(
            "gcell.epigenome.chipatlas.metadata.get_db_path", return_value=db_path
        ):
            meta = ChipAtlasMetadata.__new__(ChipAtlasMetadata)
            meta.db_path = db_path
            meta.max_workers = 4

            assert meta._get_stored_mode() == "lite"
            assert meta.is_lite_mode() is True

    def test_chipatlas_metadata_mode_property(self, tmp_path):
        """Test ChipAtlas metadata_mode and is_lite_mode properties."""
        db_path = tmp_path / "test_chipatlas.db"

        # Create mock database
        with sqlite3.connect(db_path) as conn:
            conn.executescript(
                """
                CREATE TABLE experiments (experiment_id TEXT PRIMARY KEY, assembly TEXT);
                CREATE TABLE files (file_name TEXT PRIMARY KEY, assembly TEXT);
                CREATE TABLE antigens (id INTEGER PRIMARY KEY, assembly TEXT);
                CREATE TABLE celltypes (id INTEGER PRIMARY KEY, assembly TEXT);
                CREATE TABLE metadata_info (key TEXT PRIMARY KEY, value TEXT);
            """
            )
            conn.execute(
                "INSERT INTO metadata_info VALUES (?, ?)",
                ("metadata_mode", "full"),
            )

        with (
            patch(
                "gcell.epigenome.chipatlas.metadata.get_db_path", return_value=db_path
            ),
            patch.object(ChipAtlasMetadata, "_db_exists", return_value=True),
        ):
            ca = ChipAtlas(cache_dir=tmp_path)
            ca.metadata.db_path = db_path

            assert ca.metadata_mode == "full"
            assert ca.is_lite_mode is False
