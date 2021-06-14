#include <iostream>
#include <map>
#include <omp.h>
#include <vector>

#include "image.h"
#include "stats.h"

class ImageGroup {

public:

	void run();
	void readPairs( char *fileName ); // read pair file from match

	bool printStats; // print stats at each iteration
	bool printLinear; // print linear transforms at each iteration

	int linearIterations; // number of iterations for linear registration
	bool useScale; // use scale for linear registration
	int deformableLevels; // number of deformable levels
	int deformableIterations; // number of iteration for each level
	float linearAlpha; // update ratio
	float deformableAlpha; // update ratio
	int statIntervalUpdate; // update stats every n iterations
	float initialGridSize;
	float boundingBoxMargin;
	float inlierThreshold; // to filter outliers completely
	bool guaranteeDiffeomorphism; // to guarantee that transforms will be diffeomorpic
	float maxDisplacementRatio; // max ratio to guarantee diffeomorphism
	bool invertLandmarksCoordinates; // if true the landmarks coordinates x and y will be inverted when loading

	const char* outputFileName = "measures.csv"; // output filename of measure.csv
	void writeLinksDistances(); // write distances and probabilities between pairs to file
	bool writePairs; // if true the pairs and their distances will be written to pairs.csv

	int numberOfFixedImages; // number of already registered images
	bool useRANSAC; // use RANSAC when registering with fixed images
	int numberOfRANSACIterations; // number of RANSAC iterations
	float RANSACInlierDistance; // maximum inlier distance for RANSAC
	char *fixedTransformsDirectory;

	void readLandmarks( const char *path );

	ImageGroup() {
			boundingBoxMargin = 0.1;
			deformableAlpha = 0.02;
			deformableIterations = 200;
			fixedTransformsDirectory = 0;
			guaranteeDiffeomorphism = true;
			invertLandmarksCoordinates = true;
			initialGridSize = 100;
			inlierThreshold = 0.5;
			linearAlpha =  0.5;
			linearIterations = 50;
			maxDisplacementRatio = 0.4;
			deformableLevels = 3;
			numberOfFixedImages = 0;
			numberOfRANSACIterations = 5000;
			printLinear = false;
			printStats = false;
			RANSACInlierDistance = 50;
			statIntervalUpdate = 10;
			useRANSAC = true;
			useScale = true;
			writePairs = false;
		};

protected:

	class Measure {

	public:

		float E, landmarkAv, landmarkMax, landmarkSTD;
		Measure() : landmarkAv( 0 ), landmarkMax( 0 ), landmarkSTD( 0 ),
			E ( 0 ) {};

	};


	struct Landmark {

		int image;
		float xyz[ 3 ];

	};
	typedef std::vector < Landmark > Landmarks;

	std::vector < Image > images;
	std::vector < Measure > measures;

	void setupStats();

	void setupDeformableTransforms( int level );
	void setupLinearTransforms();

	typedef pair< int, vtkMatrix4x4* > RANSACResult;
	void RANSAC( int image );
	RANSACResult RANSACBatch( int image, int nIterations, int batch );

	void transformPoints( bool apply = false );

	void updateStats(); // compute Maxwell distribution parameters
	void countInliers(); // displays number of inliers/outliers

	double updateDeformableTransforms(); // returns error value or -1 if diffeomorphism is not guaranteed
	double updateLinearTransforms();

	void displayStats();
	void displayLinearTransforms();

	void saveIndividualDistanceHistograms();
	void saveDistanceHistograms( const char *file );
	void saveMeasures( const char *file );
	void getBoundingBox( float *box, bool all = false );

	void saveTransforms();
	void saveBoundingBox();

	// reference landmarks
	std::map < std::string, Landmarks* > landmarks;
	bool computeLandmarkDistances( Measure &measure );
	void saveLandmarkDistances();


	void readAndApplyFixedImagesTransforms();

};
