#!/bin/bash
#
# Run to reproduce the full pipeline `./exploration_pipeline.sh`

# Step 1 - Align compound metadata between the two assays
cd Profiles_level4
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute align_MOA_L1000_CellPainting.ipynb

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb

# Step 2 - Analyzing level 4 profiles

# Cell Painting
cd cell_painting
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.cellpainting_calculate_cpd_median_score_spherized.ipynb

jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 2.cellpainting_calculate_SS_MAS.ipynb

jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 3.cellpainting_calculate_null_p_values.ipynb

jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 4.cellpainting_visualization.ipynb

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb

# L1000
cd ../L1000
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.L10000_calculate_cpd_median_score.ipynb

jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 2.L10000_calculate_SS_TAS.ipynb

jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 3.L10000_calculate_null_p_values.ipynb

jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 4.L10000_visualization.ipynb

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb

# L1000 and Cell Painting comparison
cd ../L1000_cellpainting_comparison

jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute L1000_CP_comparison_visualization.ipynb

jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb

# Step 3 - Analyzing level 5 profiles (Consensus signatures)

# Cell Painting
cd ../../Consensus/cell_painting
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
