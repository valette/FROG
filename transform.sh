#!/bin/bash

start=`date +%s`

function launch {

	echo "*************************************************************"
	echo "Executing : $1 "
	$1
	echo "*************************************************************"

}

# Inspired from NiftyReg script : https://sourceforge.net/p/niftyreg/git/ci/master/tree/reg-apps/groupwise_niftyreg_run.sh

if [ $# -lt 2 ]
then
	echo ""
	echo "*******************************************************************************"
	echo "Two arguments are expected to run this script:"
	echo "- File with contains parameters"
	echo "- output spacing"
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
echo ">>> There are ${IMG_NUMBER} input images to transform :"

for (( CUR_IT=0; CUR_IT<${IMG_NUMBER}; CUR_IT++ ))
do
echo ${IMG_INPUT[CUR_IT]}
done
echo "******************************************************"

# SETUP EXECUTABLES
ROOT_DIR=$(cd `dirname $0` && pwd)
AVERAGEVOLUMES=$ROOT_DIR/bin/AverageVolumes
DUMMYVOLUMEGENERATOR=$ROOT_DIR/bin/DummyVolumeGenerator
VOLUMETRANSFORM=$ROOT_DIR/bin/VolumeTransform

cd $RES_FOLDER

launch "$DUMMYVOLUMEGENERATOR bbox.json $2"
files="";

for (( CUR_IT=0; CUR_IT<${IMG_NUMBER}; CUR_IT++ ))
do
	IMG=${IMG_INPUT[CUR_IT]};
	TRANS=transformed$CUR_IT.nii.gz
	if [ ! -e "$TRANS" ]; then
		launch "$VOLUMETRANSFORM $IMG dummy.mhd -t transforms/$CUR_IT.json -o $TRANS"
	else
		echo "Transformed file $TRANS already exists, skip transformation"
	fi
	files+=" $TRANS"
done

launch "$AVERAGEVOLUMES $files"

end=`date +%s`
runtime=$((end-start))
echo "Total processing time : $runtime seconds"