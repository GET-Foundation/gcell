"""Tests for celltype module."""

import pytest


@pytest.fixture(scope="module")
def demo_loader():
    """Initialize GETDemoLoader once for all tests in this module."""
    from gcell.cell.celltype import GETDemoLoader

    return GETDemoLoader()


def test_get_demo_loader_init(demo_loader):
    """Test that GETDemoLoader initializes correctly."""
    assert demo_loader is not None
    assert hasattr(demo_loader, "available_celltypes")
    assert hasattr(demo_loader, "cfg")
    assert hasattr(demo_loader, "s3_file_sys")


def test_get_demo_loader_available_celltypes(demo_loader):
    """Test that available celltypes are loaded."""
    celltypes = demo_loader.available_celltypes
    assert isinstance(celltypes, list)
    assert len(celltypes) > 0
    # Check some expected celltypes exist
    assert "Plasma Cell" in celltypes or len(celltypes) > 10


def test_get_demo_loader_cell_type_mappings(demo_loader):
    """Test that cell type ID/name mappings are available."""
    assert hasattr(demo_loader, "cell_type_id_to_name")
    assert hasattr(demo_loader, "cell_type_name_to_id")
    assert len(demo_loader.cell_type_id_to_name) > 0
    assert len(demo_loader.cell_type_name_to_id) > 0


@pytest.mark.slow
def test_load_celltype(demo_loader):
    """Test loading a celltype from S3.

    This test is marked as slow because it downloads data from S3.
    """
    # Pick the first available celltype
    celltype_name = demo_loader.available_celltypes[0]
    cell = demo_loader.load_celltype(celltype_name)

    assert cell is not None
    assert hasattr(cell, "gene_annot")
    assert hasattr(cell, "preds")
    assert hasattr(cell, "obs")
    assert len(cell.gene_annot) > 0
