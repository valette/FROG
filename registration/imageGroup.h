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

struct Link {

	unsigned short image1;
	unsigned short point1;

	unsigned short image2;
	unsigned short point2;

	float difference[ 3 ]; // point1 - point2
	float distance; // | point1 - point2 |

	float weight; // inlier probability

};


class Measure {

public:

	float E, landmarkAv, landmarkMax, landmarkSTD;
	Measure() : landmarkAv( 0 ), landmarkMax( 0 ), landmarkSTD( 0 ),
		E ( 0 ) {};

};

class ImageGroup {

public:
	
	std::vector < Image > images;
	std::vector< Link > links; 

	std::vector < Measure > measures;

	void run();

	void setupStats();

	void setupDeformableTransforms( int level );
	void setupLinearTransforms();

	void transformPoints( bool apply = false );
	void updateDistances(); // compute all pair distances

	void updateStats(); // compute Maxwell distribution parameters

	float updateWeights(); // compute inlier probabilities, returns average distance

	bool updateDeformableTransforms(); // return true if diffeomorphism is guaranteed
	void updateLinearTransforms();
	void averageLinearTransforms1(); // set volume 0 as reference
	void averageLinearTransforms2(); // compute average scale

	void displayStats();
	void displayLinearTransforms();

	void saveIndividualDistanceHistograms();
	void saveDistanceHistograms( const char *file );
	void saveMeasures( const char *file );

	bool printStats; // print stats at each iteration
	bool printLinear; // print linear transforms at each iteration

	void readPairs( char *fileName ); // read pair file from match
	void getBoundingBox( float *box );

	int linearIterations; // number of iterations for linear registration
	bool useScale; // use scale for linear registration
	int deformableLevels; // number of deformable levels
	int deformableIterations; // number of iteration for each level
	float linearAlpha; // update ratio
	int linearAverageMethod; // 0 : set volume 0 scale to 1. 1 : use scale average
	float deformableAlpha; // update ratio
	int statIntervalUpdate; // update stats every n iterations
	float initialGridSize;
	float boundingBoxMargin;
	float inlierThreshold; // to filter outliers completely
	bool guaranteeDiffeomorphism; // to guarantee that transforms will be diffeomorpic
	float maxDisplacementRatio; // max ratio to guarantee diffeomorphism

	void saveTransforms();
	void saveBoundingBox();

	// reference landmarks
	std::map < std::string, Landmarks* > landmarks;
	void readLandmarks( const char *path );
	bool computeLandmarkDistances( Measure &measure );

	ImageGroup() : linearIterations( 200 ), deformableIterations( 200 ),
		linearAlpha( 0.5 ), linearAverageMethod( 1 ), deformableAlpha( 0.02 ),
		deformableLevels( 3 ), guaranteeDiffeomorphism( true ),
		maxDisplacementRatio( 0.4 ), useScale( true ), statIntervalUpdate( 10 ),
		initialGridSize( 100 ), boundingBoxMargin( 0.1 ), inlierThreshold( 0.5 ),
		printStats( false ), printLinear( false ) {};

};
