#!/bin/bash

start=`date +%s`

function launch {

	echo "*************************************************************"
	echo "Executing : $1 "
	$1
	RET=$?
	echo "*************************************************************"
	return $RET
}

# Inspired from NiftyReg script : https://sourceforge.net/p/niftyreg/git/ci/master/tree/reg-apps/groupwise_niftyreg_run.sh

if [ $# -lt 1 ]
then
	echo ""
	echo "*******************************************************************************"
	echo "At least one argument is expected to run this script:"
	echo "- File with contains parameters"
	echo "- spacing (optional)"
	echo "example: $0 params.sh 2 "
	echo "*******************************************************************************"
	echo ""
	exit
fi

# read input parameters
. $1

# check arguments

if [ ${#IMG_INPUT[@]} -lt 2 ]
then
	echo "Less than 2 images have been specified"
	echo "Exit ..."
	exit
fi

IMG_NUMBER=${#IMG_INPUT[@]}

echo ""
echo "******************************************************"
echo ">>> There are ${IMG_NUMBER} input images :"

for (( CUR_IT=0; CUR_IT<${IMG_NUMBER}; CUR_IT++ ))
do
echo ${IMG_INPUT[CUR_IT]}
done
echo "******************************************************"

# SETUP EXECUTABLES
ROOT_DIR=$(cd `dirname $0` && pwd)
CHECKDIFFEOMORPHISM=$ROOT_DIR/bin/CheckDiffeomorphism

cd $RES_FOLDER
N=0;
SPACING=${2:--1}

for (( CUR_IT=0; CUR_IT<${IMG_NUMBER}; CUR_IT++ ))
do
	IMG=${IMG_INPUT[CUR_IT]};
	TRANS=transformed$CUR_IT.nii.gz
	launch "$CHECKDIFFEOMORPHISM $IMG tfm/$CUR_IT.tfm $SPACING"
	if [ $? -ne 0 ]
	then
		N=1
	fi
done

end=`date +%s`
runtime=$((end-start))
echo "Total processing time : $runtime seconds"

if [ $N -eq 0 ]
then
	echo "All transforms are diffeomorphic"
else
	echo "Some transforms are not diffeomorphic. Please check log"
fi
