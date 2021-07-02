# LINCS profiling complimentarity

Assessing information contained in different profiling data modalities.

The two assays we evaluated are:

* L1000 (gene expression)
* Cell Painting (morphology)

## Evaluation

First, we aligned profiles with [`1.Data-exploration/Profiles_level4/align_MOA_L1000_CellPainting.ipynb`](1.Data-exploration/Profiles_level4/align_MOA_L1000_CellPainting.ipynb).
This resulted in us using only the perturbations measured in common (by Broad ID and dose).
Next, using these perturbations, we explored various data levels and evaluated the assays in the following ways:

- Level 4 spherized feature selected profiles
  - Median pairwise replicate reproducibility compared to a carefully matched null distribution
  - Morphological activity score and signature strength
  - PCA and sample clustering
  - Multi-label / Multi-class mechanism of action (MOA) prediction
    - We used the top models as presented in a recent [Kaggle competition](https://www.kaggle.com/c/lish-moa)
- Level 5 consensus profiles
  - MOA median pairwise similarity compared to a null distribution
    - Across increasing treatment dose
  - MOA query recall

## Computational environment

We use a combination of conda and pip to manage the proper python packages for data assessment and model predictions.
To reproduce our environment run the following:

```bash
# In the top folder of the directory
conda env create --force --file environment.yml && conda activate lincs-complimentarity && cd 2.MOA-prediction/ && python setup.py && cd ..
```

We also need to setup custom computational environments for tensorflow and pytorch for the MOA prediction analysis.

```bash
# Navigate into the MOA prediction folder
cd 2.MOA-prediction

# Step 1 - Tensorflow
# Initialize a virtual environment
python3 -m venv tensorflow_env

# Activate the environment
source tensorflow_env/bin/activate

# Upgrade pip if necessary
# python3 -m pip install --upgrade pip

# Install tensorflow requirements
python3 -m pip install -r tensorflow_requirements.txt

# Step 2 - Pytorch
python3 -m venv pytorch_env
source pytorch_env/bin/activate
python3 -m pip install -r pytorch_requirements.txt && python3 setup.py
```
