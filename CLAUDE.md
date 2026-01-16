# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

gcell is a Python library for genomic data analysis, focusing on cell type-specific regulatory analysis, DNA sequence manipulation, protein structure prediction, and pathway analysis. It serves as a toolkit for AI agents to discover new biology through predictive models.

## Build & Development Commands

```bash
# Environment setup
uv sync                        # Install dependencies
uv sync --all-extras           # Install with all optional deps
pre-commit install             # Install git hooks

# Testing
uv run pytest                  # Run all tests
uv run pytest -n auto          # Run tests in parallel
uv run pytest -k "pattern"     # Run specific tests by pattern
uv run pytest tests/gcell/dna/ # Run tests in specific directory

# Code quality
ruff check --fix .             # Fix linting issues
ruff format .                  # Format code
pre-commit run --all-files     # Run all pre-commit hooks

# Documentation
uv run sphinx-build docs docs/_build  # Build Sphinx docs

# Running scripts
uv run python script.py        # Run Python script in env
```

## Architecture

```
src/gcell/
├── cell/       # Cell type analysis (GETCellType, GETHydraCellType)
│               # Gene expression, motif analysis, causal networks
├── dna/        # DNA sequences (Genome, Sequence, Motif, Track)
│               # Motif scanning, peak files, genomic regions
├── rna/        # Gene annotations (Gencode, GTF, Gene, TSS)
├── protein/    # Protein analysis (AF2, UniProt, InterPro, STRING)
├── ontology/   # Pathway analysis via gprofiler (GO, KEGG, Reactome)
├── utils/      # Causal inference (lingam), PDB viewer, S3 support
└── config/     # Hydra/OmegaConf configuration files
```

### Data Management
- Data stored in `~/.gcell_data/` with subdirs: `annotations/`, `genomes/`, `cache/`
- Environment variables: `GCELL_*` for directory overrides
- Uses `pooch` for cached downloads, `s3fs` for cloud storage

### Key Dependencies
- **torch is optional**: Install with `pip install gcell[torch]`
- **zarr < 3.0**: Required due to S3 compatibility issues

## Code Style

- Line length: 88 characters
- Formatter/Linter: Ruff
- Docstrings: NumPy format
- Imports: isort with `gcell` as first-party

## Testing

- Framework: pytest with doctest support
- Use module-scoped fixtures for expensive operations (data loading)
- Image comparison tests for plotting: store expected images in `tests/_images/`

## Loading Pre-inferred Cell Types

```python
from gcell.cell.celltype import GETDemoLoader
g = GETDemoLoader()
print(g.available_celltypes)
g.load_celltype('Plasma Cell')
```
