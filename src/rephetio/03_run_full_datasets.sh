#!/bin/bash

# see .02 for parallel edits made

LOAD_DIR='../out/' # line adjusted from original

for folder in ${LOAD_DIR}20*
    do python3 full_dataset_training.py $folder
done
