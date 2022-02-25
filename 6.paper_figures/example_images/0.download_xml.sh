#!/bin/bash

# Step 0: Setup AWS CLI

# Step 1: Download XML files
aws s3 cp --recursive --no-sign-request \
  s3://cellpainting-gallery/lincs/broad/images/2016_04_01_a549_48hr_batch1/images/ xml/ \
  --exclude "*" --include "*.xml"
