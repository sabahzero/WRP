#!/bin/bash

LOAD_DIR='out/' # line adjusted from original

for folder in ${LOAD_DIR}* # line adjusted from original (adjusted to output folders)
    do python3 n_fold_CV_training.py $folder # do needs to be on the same line when in terminal, need to make sure it's python3 (edit)
done

# need to be able to debug code
# mike had 2018-04, 2018-02, etc... with nodes and edges .csvs within these (that's why there was a 20*)

# When attempted to run, continued to get
##  Traceback (most recent call last):
##   File "n_fold_CV_training.py", line 11, in <module>
##     from glmnet import LogitNet # Need to install pip install glmnet
##   File "/Users/sulhasan/opt/anaconda3/lib/python3.8/site-packages/glmnet/__init__.py", line 3, in <module>
##     from .logistic import LogitNet
##   File "/Users/sulhasan/opt/anaconda3/lib/python3.8/site-packages/glmnet/logistic.py", line 14, in <module>
##     from _glmnet import lognet, splognet, lsolns
## ImportError: dlopen(/Users/sulhasan/opt/anaconda3/lib/python3.8/site-packages/_glmnet.cpython-38-darwin.so, 2): Library not loaded: /usr/local/opt/gcc/lib/gcc/9/libgfortran.5.dylib
##   Referenced from: /Users/sulhasan/opt/anaconda3/lib/python3.8/site-packages/_glmnet.cpython-38-darwin.so
##   Reason: image not found

## Tried, and didn't work
### sudo pip install --upgrade --force-reinstall glmnet
### conda install python$pythonversion$
### conda install -c conda-forge glmnet
### conda install -c conda-forge/label/gcc7 glmnet
### conda install -c conda-forge/label/cf201901 glmnet
### conda install -c conda-forge/label/cf202003 glmnet

### git clone git@github.com:civisanalytics/python-glmnet.git
### "Exception: Failed to find libgfortran.3.dylib"
#### Installed Fortran, that didn't work
#### https://github.com/eth-cscs/abcpy/issues/43

###### Might? have worked after doing this (install fortran), error went away
# http://hpc.sourceforge.net/

# ModuleNotFoundError: No module named 'glmnet'
# import glmnet_python added to respective scripts but doesnt do anything
# https://github.com/civisanalytics/python-glmnet doesn't do anything