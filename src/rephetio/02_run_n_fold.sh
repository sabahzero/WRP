#!/bin/bash

LOAD_DIR='../out/' # line adjusted from original

for folder in ${LOAD_DIR}20*
    do python n_fold_CV_training.py $folder
done