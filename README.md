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
├── data/                   # Data files and analysis results
│   ├── AlignedModels/      # ✅ Aligned anatomical models (CSV format) - 88 files
│   ├── ArotaModels/        # 🔍 Original STL models (7-15MB each) - 88 files
│   ├── Final/              # 🔍 Generated PLY point clouds - anatomical_*/
│   │   ├── anatomical_*/   # Individual model directories
│   │   │   ├── *_original.ply
│   │   │   ├── *_PC1_minus2sigma.ply
│   │   │   ├── *_PC1_plus2sigma.ply
│   │   │   └── ... (PC2, PC3 variants)
│   └── Evaluate_output/    # ✅ Analysis results and evaluations
│       ├── evaluation_generalization_loormse.csv
│       ├── PCA_eigenvalues_analysis.txt
│       └── specificity_out/    # K1-K10 specificity analysis
└── README.md               # This file

✅ = Included in repository
🔍 = Large files, contact for dataset access
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

## Dataset Information

### Included Data ✅
- **88 Aligned Models**: CSV format anatomical models ready for PCA analysis
- **Evaluation Results**: Complete analysis outputs including RMSE, specificity metrics
- **PCA Analysis**: Eigenvalue analysis and component summaries

### Large Dataset Access 🔍
The complete dataset includes:
- **88 STL Models**: Original anatomical aorta models (7-15MB each, ~800MB total)
- **616 PLY Files**: Generated point clouds with PC variations (±2σ for PC1-3)
- **Additional Outputs**: Generated meshes and intermediate processing files

*For access to the complete dataset, please contact the repository owner.*

## Notes

- Repository optimized for code sharing and reproducibility
- Large binary files managed separately to maintain reasonable clone size
- All analysis can be reproduced with the provided aligned CSV data
- Build artifacts and IDE files are excluded via .gitignore

## Academic Context

This project was developed for academic research in statistical shape modeling of anatomical structures, specifically focusing on aorta analysis using advanced computational geometry techniques.
