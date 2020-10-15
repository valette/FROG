#include <iostream>
#include <map>
#include <vector>

#include "image.h"
#include "stats.h"

struct Landmark {

	int image;
	float xyz[ 3 ];

};

typedef std::vector < Landmark > Landmarks;

class Measure {

public:

	float E, landmarkAv, landmarkMax, landmarkSTD;
	Measure() : landmarkAv( 0 ), landmarkMax( 0 ), landmarkSTD( 0 ),
		E ( 0 ) {};

};

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

	int numberOfFixedImages;
	char *fixedTransformsDirectory;

	void readLandmarks( const char *path );

	ImageGroup() : linearIterations( 50 ), deformableIterations( 200 ),
		linearAlpha( 0.5 ), deformableAlpha( 0.02 ),
		deformableLevels( 3 ), fixedTransformsDirectory( 0 ),guaranteeDiffeomorphism( true ),
		invertLandmarksCoordinates( true ), maxDisplacementRatio( 0.4 ),
		numberOfFixedImages( 0 ), useScale( true ), statIntervalUpdate( 10 ),
		initialGridSize( 100 ), boundingBoxMargin( 0.1 ), inlierThreshold( 0.5 ),
		printStats( false ), printLinear( false ) {};

protected:

	std::vector < Image > images;
	std::vector < Measure > measures;

	void setupStats();

	void setupDeformableTransforms( int level );
	void setupLinearTransforms();

	void transformPoints( bool apply = false );

	void updateStats(); // compute Maxwell distribution parameters

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
