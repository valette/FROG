#!/bin/sh

start=`date +%s`

function launch {

	echo "*************************************************************"
	echo "Executing : $1 "
	$1
	echo "*************************************************************"

}

# Inspired from NiftyReg script : https://sourceforge.net/p/niftyreg/git/ci/master/tree/reg-apps/groupwise_niftyreg_run.sh

if [ $# -lt 1 ]
then
	echo ""
	echo "*******************************************************************************"
	echo "One argument is expected to run this script:"
	echo "- File with contains parameters"
	echo "example: $0 params.sh "
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
echo ">>> There are ${IMG_NUMBER} input images to register :"

for (( CUR_IT=0; CUR_IT<${IMG_NUMBER}; CUR_IT++ ))
do
echo ${IMG_INPUT[CUR_IT]}
done
echo "******************************************************"

# SETUP EXECUTABLES
ROOT_DIR=$(cd `dirname $0` && pwd)
SURF=$ROOT_DIR/bin/surf3d
MATCH=$ROOT_DIR/bin/match
REG=$ROOT_DIR/bin/frog

# CREATE THE RESULT FOLDER
if [ ! -d ${RES_FOLDER} ]
then
	echo "The output image folder (${RES_FOLDER}) does not exist"
	mkdir ${RES_FOLDER}
	if [ ! -d ${RES_FOLDER} ]
	then
		echo "Unable to create the ${RES_FOLDER} folder"
		echo "Exit ..."
		exit
	else
		echo "The output image folder (${RES_FOLDER}) has been created"
	fi
fi

cd $RES_FOLDER
CSVFile=points.csv;

if [ -f ${CSVFile} ]
then
    rm ${CSVFile}
fi

for (( CUR_IT=0; CUR_IT<${IMG_NUMBER}; CUR_IT++ ))
do
	IMG=${IMG_INPUT[CUR_IT]};
	OUTPUT_POINTS=$RES_FOLDER/points$CUR_IT
	launch "$SURF $IMG -t 0 -n 20000 -o $OUTPUT_POINTS -s $SPACING -t $THRESHOLD -n $NPOINTS"
	echo ${OUTPUT_POINTS}.csv.gz,0,0,0 >> ${CSVFile}

done

OUTPUT_PAIRS=pairs.bin
launch "$MATCH $CSVFile -o $OUTPUT_PAIRS -d $MAX_DISTANCE"

if [ -n "$LANDMARKS" ]
then
	LOPTIONS="-l ${LANDMARKS}"
else
	LOPTIONS=""
fi

launch "$REG $OUTPUT_PAIRS $LOPTIONS"

end=`date +%s`
runtime=$((end-start))
echo "Total registration time : $runtime seconds"