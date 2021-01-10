#!/bin/sh

# Inspired by Niftireg script : https://sourceforge.net/p/niftyreg/git/ci/master/tree/reg-apps/groupwise_niftyreg_params.sh

# Array that contains the input images to create the atlas
export IMG_INPUT=(`ls ~/path/to/images/*.nii.gz`)

# folder where the results will be written
export RES_FOLDER=`pwd`/output

# SURF parameters
export SPACING=0.75
export THRESHOLD=0;
export NPOINTS=20000;
export SURF_OTHER_PARAMS="" # other possible parameters, see documentation here: https://github.com/valette/vtkOpenSURF3D

# MATCH parameters
export MAX_DISTANCE=1
export MATCH_OTHER_PARAMS="" # other possible match parameters

# FROG parameters

#export LANDMARKS=/path/to/landmarks # path to landmarks to measure registration quality
export REGISTRATION_OTHER_PARAMS=""
