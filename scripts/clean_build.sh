#!/bin/bash
# Clean build artifacts and temporary files

set -e

echo "Cleaning build artifacts..."

# Remove build directory
if [ -d "build" ]; then
    echo "  Removing build/ directory"
    rm -rf build/
fi

# Remove dist directory
if [ -d "dist" ]; then
    echo "  Removing dist/ directory"
    rm -rf dist/
fi

# Remove egg-info directories
find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true

# Remove temporary setuptools files (=8, =64, etc.)
find . -maxdepth 1 -type f -name "=*" -exec rm -f {} + 2>/dev/null || true

# Remove __pycache__ directories
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
find . -type f -name "*.pyc" -delete 2>/dev/null || true
find . -type f -name "*.pyo" -delete 2>/dev/null || true

echo "✓ Cleanup complete!"

