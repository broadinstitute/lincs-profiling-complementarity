#!/bin/bash

# Step 0 - Download and process data
cd 0.download_cellpainting_L1000_data
python download_feature_select_cp_lvl4.py --data_dir "data"
python download_feature_select_L1000_lvl4.py --data_dir "data"
cd ..

# Step 1 - Split compounds prior to splitting profiles:
# This splitting will inform training and test splits
cd 1.compound_split_train_test
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute split_cpds_train_test.ipynb
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb
cd ..

# Step 2 - Split profiles into training and test sets
cd 2.data_split
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.cellpainting_L1000_data_train_test_split.ipynb

jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.merge_cp_L1000_and_train_test_split.ipynb
        
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb
cd ..

# Step 3 - Perform MOA predictions
cd 3.moa_prediction_models
./run_models.sh
cd ..

# Step 4 - Blend predictions, acquire performance metrics, and visualize results
cd 4.model_viz
jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 0.blend_test_predictions.ipynb

jupyter nbconvert --to=html \
        --FilesWriter.build_directory=scripts/html \
        --ExecutePreprocessor.kernel_name=python3 \
        --ExecutePreprocessor.timeout=10000000 \
        --execute 1.moa_predictions_visualization.ipynb
        
jupyter nbconvert --to=script --FilesWriter.build_directory=scripts/nbconverted *.ipynb
cd ..