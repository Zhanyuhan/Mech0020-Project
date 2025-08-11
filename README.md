# MyScalismo Project - Aorta Analysis with PCA

## Overview
This project implements Principal Component Analysis (PCA) for anatomical aorta models using the Scalismo framework. It provides tools for statistical shape analysis, model evaluation, and format conversion between different 3D file formats.

## Key Features

### Core Analysis Tools
- **CombinedPCA.scala**: Main PCA implementation with comprehensive shape analysis
- **Evaluate.scala**: Model evaluation using leave-one-out cross-validation and RMSE metrics
- **Evaluate_Specificity.scala**: Specificity analysis for model validation
- **NewPrepareToCSV.scala**: Data preprocessing and CSV format conversion
- **PCA_PCratio.scala**: Principal component ratio analysis

### File Format Conversion
- **plytostlsmooth.py**: Advanced PLY to STL conversion with:
  - Rigid ICP alignment using SVD-Umeyama algorithm
  - Laplacian mesh smoothing with boundary preservation
  - Non-destructive processing (preserves original files)

## Project Structure

```
MyScalismoProject/
├── *.scala                 # Scala source files for PCA analysis
├── *.py                    # Python utilities for mesh processing
├── build.sbt               # SBT build configuration
├── data/                   # Data files (excluded from git due to size)
│   ├── AlignedModels/      # Aligned anatomical models (CSV format)
│   ├── ArotaModels/        # Original STL models
│   ├── Final/              # Generated PLY point clouds
│   └── Evaluate_output/    # Analysis results and evaluations
└── README.md               # This file
```

## Technologies Used

- **Scalismo**: Swiss Statistical Shape Modeling framework
- **Scala 3**: Modern functional programming language
- **Python**: For mesh processing and format conversion
- **NumPy/SciPy**: Scientific computing libraries
- **SBT**: Scala Build Tool

## Key Algorithms

### PCA Analysis (Scala)
- Statistical shape modeling of anatomical structures
- Eigenvalue decomposition for principal components
- Cross-validation with leave-one-out methodology
- Specificity and generalization analysis

### Mesh Processing (Python)
- **Iterative Closest Point (ICP)**: Rigid alignment with scale
- **Umeyama Algorithm**: SVD-based registration
- **Laplacian Smoothing**: Umbrella operator with boundary preservation
- **Format Conversion**: PLY ↔ STL with topology preservation

## Data Pipeline

1. **Input**: Original STL anatomical models
2. **Preprocessing**: Alignment and normalization
3. **PCA Analysis**: Statistical shape modeling
4. **Evaluation**: Cross-validation and metrics
5. **Output**: Generated models and analysis results

## Usage

### Running PCA Analysis
```scala
// Compile and run using SBT
sbt compile
sbt "runMain CombinedPCA"
```

### Mesh Processing
```python
# PLY to STL conversion with smoothing
python plytostlsmooth.py --ply_dir "data/Final" --stl_dir "data/ArotaModels" --output_dir "output" --all
```

## Code Quality

- **Fully English**: All Chinese comments and messages have been translated to English
- **Professional Standards**: Clean, readable code following best practices
- **No Emojis**: Professional presentation suitable for academic/research use
- **Comprehensive Documentation**: Detailed function documentation and comments

## Requirements

### Scala Dependencies
- Scalismo framework
- Scala 3.x
- SBT 1.x

### Python Dependencies
- NumPy
- SciPy (optional, for optimized nearest neighbor)
- Python 3.x

## Notes

- Large data files (STL, PLY, CSV) are excluded from git due to size limitations
- Build artifacts and IDE files are also excluded via .gitignore
- Original functionality is preserved while improving code readability

## Academic Context

This project was developed for academic research in statistical shape modeling of anatomical structures, specifically focusing on aorta analysis using advanced computational geometry techniques.
