#!/bin/bash
#
# Run to reproduce the full pipeline `./clustering_pca_pipeline.sh`

# Step 0 - Perform the clustering and PCA analysis in Cell Painting data
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.cp_pca_clustering_analysis.ipynb

# Step 1 - Perform the clustering and PCA analysis in L1000 data
cd Profiles_level4
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.cp_pca_clustering_analysis.ipynb

# Step 2 - Visualize clustering metrics
cd Profiles_level4
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 2.clustering_visualization.ipynb

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb
