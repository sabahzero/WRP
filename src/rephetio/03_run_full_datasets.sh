#!/bin/bash

LOAD_DIR='../out/' # line adjusted from original

for folder in ${LOAD_DIR}20*
    do python full_dataset_training.py $folder
done
