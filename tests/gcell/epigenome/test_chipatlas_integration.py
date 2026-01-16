"""Integration tests for the ChIP-Atlas query module.

These tests use real data from ChIP-Atlas but with limited rows
to keep the tests fast and bandwidth-efficient.
"""

import pandas as pd
import pytest

from gcell.epigenome.chipatlas import ChipAtlas


@pytest.fixture(scope="module")
def chipatlas_instance(tmp_path_factory):
    """Create a ChipAtlas instance with limited data for testing.

    Uses max_rows=100 to download only ~100 rows from each metadata file,
    which is enough to test the functionality without downloading the full
    database (~500MB+).
    """
    cache_dir = tmp_path_factory.mktemp("chipatlas_test")
    # Download only 100 rows for testing (uses HTTP Range requests)
    ca = ChipAtlas(
        cache_dir=cache_dir,
        force_refresh=True,
        max_rows=100,
    )
    return ca


class TestChipAtlasIntegration:
    """Integration tests with real ChIP-Atlas data."""

    def test_metadata_downloaded(self, chipatlas_instance):
        """Test that metadata was downloaded successfully."""
        ca = chipatlas_instance
        assert ca.metadata.db_path.exists()

    def test_search_experiments(self, chipatlas_instance):
        """Test searching for experiments."""
        ca = chipatlas_instance
        # Search without filters - should return some results
        results = ca.search(limit=10)
        assert isinstance(results, pd.DataFrame)
        assert len(results) <= 10
        if len(results) > 0:
            assert "experiment_id" in results.columns
            assert "assembly" in results.columns
            assert "antigen" in results.columns

    def test_search_by_assembly(self, chipatlas_instance):
        """Test searching by assembly."""
        ca = chipatlas_instance
        # Get available assemblies from data
        assemblies = ca.get_assemblies()
        assert isinstance(assemblies, list)

        if len(assemblies) > 0:
            # Search for experiments in first available assembly
            results = ca.search(assembly=assemblies[0], limit=5)
            assert isinstance(results, pd.DataFrame)
            if len(results) > 0:
                assert all(results["assembly"] == assemblies[0])

    def test_get_antigens(self, chipatlas_instance):
        """Test getting available antigens."""
        ca = chipatlas_instance
        antigens = ca.get_antigens()
        assert isinstance(antigens, pd.DataFrame)
        if len(antigens) > 0:
            assert "antigen" in antigens.columns
            assert "assembly" in antigens.columns

    def test_get_celltypes(self, chipatlas_instance):
        """Test getting available cell types."""
        ca = chipatlas_instance
        celltypes = ca.get_celltypes()
        assert isinstance(celltypes, pd.DataFrame)
        if len(celltypes) > 0:
            assert "cell_type" in celltypes.columns
            assert "assembly" in celltypes.columns

    def test_get_celltype_classes(self, chipatlas_instance):
        """Test getting cell type classes."""
        ca = chipatlas_instance
        classes = ca.get_celltype_classes()
        assert isinstance(classes, list)

    def test_get_antigen_classes(self, chipatlas_instance):
        """Test getting antigen classes."""
        ca = chipatlas_instance
        classes = ca.get_antigen_classes()
        assert isinstance(classes, list)

    def test_summary(self, chipatlas_instance):
        """Test getting summary statistics."""
        ca = chipatlas_instance
        summary = ca.summary()
        assert isinstance(summary, pd.DataFrame)
        if len(summary) > 0:
            assert "assembly" in summary.columns
            assert "experiment_count" in summary.columns

    def test_search_as_experiments(self, chipatlas_instance):
        """Test search returning experiment objects."""
        ca = chipatlas_instance
        experiments = ca.search(limit=3, as_experiments=True)
        assert isinstance(experiments, list)
        if len(experiments) > 0:
            exp = experiments[0]
            assert hasattr(exp, "experiment_id")
            assert hasattr(exp, "assembly")
            assert hasattr(exp, "bigwig_url")
            # Check URL generation
            assert exp.bigwig_url.startswith("https://")
            assert exp.experiment_id in exp.bigwig_url


class TestChipAtlasPeakDownload:
    """Tests for downloading peak data (requires network)."""

    @pytest.fixture
    def real_experiment_id(self, chipatlas_instance):
        """Get a real experiment ID for testing downloads."""
        ca = chipatlas_instance
        results = ca.search(limit=1)
        if len(results) > 0:
            return results.iloc[0]["experiment_id"], results.iloc[0]["assembly"]
        pytest.skip("No experiments available in test data")

    @pytest.mark.slow
    def test_get_peaks(self, chipatlas_instance, real_experiment_id):
        """Test downloading peak data for a real experiment."""
        ca = chipatlas_instance
        exp_id, assembly = real_experiment_id

        try:
            peaks = ca.get_peaks(exp_id, assembly=assembly, threshold=5)
            assert isinstance(peaks, pd.DataFrame)
            if len(peaks) > 0:
                assert "Chromosome" in peaks.columns
                assert "Start" in peaks.columns
                assert "End" in peaks.columns
        except Exception as e:
            # Peak data may not be available for all experiments
            if "404" in str(e) or "Not Found" in str(e):
                pytest.skip(f"Peak data not available for {exp_id}")
            raise


class TestHTTPRangeDownload:
    """Tests specifically for HTTP Range request functionality."""

    def test_partial_download(self):
        """Test that partial downloads work correctly."""
        import urllib.request

        url = "https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab"

        # Download first 5KB
        req = urllib.request.Request(url)
        req.add_header("Range", "bytes=0-5000")

        with urllib.request.urlopen(req, timeout=30) as response:
            content = response.read()

        # Should have received partial content
        assert len(content) > 0
        assert len(content) <= 5001  # May be slightly more due to HTTP

        # Content should be valid TSV
        lines = content.decode("utf-8", errors="ignore").split("\n")
        assert len(lines) > 1

        # First line should have tab-separated fields
        first_line = lines[0]
        fields = first_line.split("\t")
        assert len(fields) >= 2  # At least experiment_id and assembly

    def test_max_rows_limits_data(self, tmp_path):
        """Test that max_rows parameter limits downloaded data."""
        from gcell.epigenome.chipatlas.metadata import ChipAtlasMetadata

        # Create metadata with very limited rows
        meta = ChipAtlasMetadata.__new__(ChipAtlasMetadata)
        meta.db_path = tmp_path / "test.db"
        meta.max_workers = 2
        meta.max_rows = 10

        # Download a single file with limit
        name, df = meta._download_metadata_file(
            "experimentList",
            "https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab",
            max_bytes=5000,  # ~10 rows
        )

        assert len(df) <= 20  # Should have limited rows
        assert "experiment_id" in df.columns
