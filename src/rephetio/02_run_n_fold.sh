#!/bin/bash

LOAD_DIR='out/' # line adjusted from original

for folder in ${LOAD_DIR}* # line adjusted from original (adjusted to output folders)
    do python n_fold_CV_training.py $folder # do needs to be on the same line when in terminal
done

# need to be able to debug code
# mike had 2018-04, 2018-02, etc... with nodes and edges .csvs within these (that's why there was a 20*)