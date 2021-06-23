#!/bin/bash
#
# Run to reproduce the full pipeline `./exploration_pipeline.sh`

jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.cellpainting_moa_median_scores_calculation.ipynb

jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 2.cellpainting_null_p_values_calculation.ipynb

jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 3.cellpainting_visualization.ipynb

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb

# L1000
cd ../L1000
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.L1000_moa_median_scores_calculation.ipynb

jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 2.L1000_null_p_values_calculation.ipynb

jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 3.L1000_visualization.ipynb

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb

# L1000 and Cell Painting comparison
cd ../L1000_cell_painting_comparison

jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute L1000_CP_comparison_visualization.ipynb

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb
