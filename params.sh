#!/bin/sh

# Inspired from Niftireg script : https://sourceforge.net/p/niftyreg/git/ci/master/tree/reg-apps/groupwise_niftyreg_params.sh

# Array that contains the input images to create the atlas
export IMG_INPUT=(`ls ~/path/to/images/*.nii.gz`)

# folder where the results will be written
export RES_FOLDER=`pwd`/registration

# SURF parameters
export SPACING=0.75
export THRESHOLD=0;
export NPOINTS=20000;

# MATCH parameters
export MAX_DISTANCE=1

# Optional path to landmarks
#export LANDMARKS=/path/to/landmarks
