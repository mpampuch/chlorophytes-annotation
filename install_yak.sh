#!/bin/bash

# Installation script for YAK and pipeline dependencies
# Usage: bash install_yak.sh [install_dir]

set -euo pipefail

# Default installation directory
INSTALL_DIR="${1:-$HOME/tools}"
YAK_DIR="$INSTALL_DIR/yak"

echo "=========================================="
echo "Installing YAK and Pipeline Dependencies"
echo "=========================================="

# Create installation directory
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

# Check if git is available
if ! command -v git &> /dev/null; then
    echo "ERROR: git is required but not installed"
    exit 1
fi

# Check if make is available
if ! command -v make &> /dev/null; then
    echo "ERROR: make is required but not installed"
    exit 1
fi

# Check if gcc/clang is available
if ! command -v gcc &> /dev/null && ! command -v clang &> /dev/null; then
    echo "ERROR: C compiler (gcc or clang) is required but not installed"
    exit 1
fi

echo "Installing YAK..."

# Clone YAK repository
if [ -d "$YAK_DIR" ]; then
    echo "YAK directory already exists, updating..."
    cd "$YAK_DIR"
    git pull
else
    echo "Cloning YAK repository..."
    git clone https://github.com/lh3/yak.git "$YAK_DIR"
    cd "$YAK_DIR"
fi

# Compile YAK
echo "Compiling YAK..."
make clean || true
make

# Check if compilation was successful
if [ -x "./yak" ]; then
    echo "YAK compiled successfully!"
    
    # Test YAK installation
    echo "Testing YAK installation..."
    ./yak 2>&1 | head -3 || echo "YAK executable found"
    
    # Add to PATH instructions
    echo ""
    echo "=========================================="
    echo "Installation Complete!"
    echo "=========================================="
    echo "YAK has been installed to: $YAK_DIR"
    echo ""
    echo "To use YAK, add it to your PATH:"
    echo "export PATH=\"$YAK_DIR:\$PATH\""
    echo ""
    echo "Add this line to your ~/.bashrc or ~/.zshrc to make it permanent:"
    echo "echo 'export PATH=\"$YAK_DIR:\$PATH\"' >> ~/.bashrc"
    echo ""
    
    # Check Nextflow installation
    echo "Checking Nextflow installation..."
    if command -v nextflow &> /dev/null; then
        NEXTFLOW_VERSION=$(nextflow -version | head -1 | awk '{print $3}')
        echo "Nextflow found: version $NEXTFLOW_VERSION"
        
        # Check if version is sufficient
        if [[ $(echo "$NEXTFLOW_VERSION 22.10.0" | tr " " "\n" | sort -V | head -n1) == "22.10.0" ]]; then
            echo "Nextflow version is sufficient (≥22.10.0)"
        else
            echo "WARNING: Nextflow version may be too old (recommend ≥22.10.0)"
        fi
    else
        echo "Nextflow not found. Please install Nextflow:"
        echo "curl -s https://get.nextflow.io | bash"
        echo "sudo mv nextflow /usr/local/bin/"
    fi
    
    # Check Python installation
    echo ""
    echo "Checking Python installation..."
    if command -v python3 &> /dev/null; then
        PYTHON_VERSION=$(python3 --version | awk '{print $2}')
        echo "Python3 found: version $PYTHON_VERSION"
    else
        echo "WARNING: Python3 not found. Required for report generation."
    fi
    
    echo ""
    echo "=========================================="
    echo "Next Steps:"
    echo "=========================================="
    echo "1. Add YAK to your PATH (see instructions above)"
    echo "2. Test the pipeline:"
    echo "   nextflow run main.nf --help"
    echo "3. Run with your data:"
    echo "   nextflow run main.nf \\"
    echo "     --illumina_reads 'data/illumina/*_{1,2}.fastq.gz' \\"
    echo "     --nanopore_reads 'data/nanopore/*.fastq.gz'"
    echo ""
    
else
    echo "ERROR: YAK compilation failed!"
    echo "Please check the error messages above."
    exit 1
fi

echo "Installation script completed successfully!"