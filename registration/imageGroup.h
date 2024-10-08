#include <iostream>
#include <map>
#include <omp.h>
#include <vector>

#include "image.h"
#include "stats.h"
#include "../vtkOpenSURF3D/picojson.h"

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
	float linearInitializationAnchor[ 3 ];
	float linearAlpha; // update ratio
	float deformableAlpha; // update ratio
	int statIntervalUpdate; // update stats every n iterations
	float initialGridSize;
	float boundingBoxMargin;
	float inlierThreshold; // to filter outliers completely
	bool guaranteeDiffeomorphism; // to guarantee that transforms will be diffeomorpic
	float maxDisplacementRatio; // max ratio to guarantee diffeomorphism
	bool invertLandmarksCoordinates; // if true the landmarks coordinates x and y will be inverted when loading
	float landmarksConstraintsWeight; // weight for landmark constraints

	const char* outputFileName = "measures.csv"; // output filename of measure.csv
	void writeLinksDistances(); // write distances and probabilities between pairs to file
	bool writePairs; // if true the pairs and their distances will be written to pairs.csv

	int numberOfFixedImages; // number of already registered images
	bool useRANSAC; // use RANSAC when registering with fixed images
	int numberOfRANSACIterations; // number of RANSAC iterations
	float RANSACInlierDistance; // maximum inlier distance for RANSAC
	float RANSACMaxScale; // maximum allowed scale for RANSAC
	char *fixedTransformsDirectory;
	bool writeSingleFileTransforms; // outputs a single big JSON file for each transform

	std::string transformSubdirectory;
	std::string errorMapsSubdirectory;
	void addLandmarks( const char *path, bool asConstraints = false );

	ImageGroup() {
			boundingBoxMargin = 0.1;
			deformableAlpha = 0.02;
			deformableIterations = 200;
			fixedTransformsDirectory = 0;
			guaranteeDiffeomorphism = true;
			invertLandmarksCoordinates = true;
			landmarksConstraintsWeight = 50;
			initialGridSize = 100;
			inlierThreshold = 0.5;
			linearAlpha = 0.5;
			linearInitializationAnchor[ 0 ] = 0.5;
			linearInitializationAnchor[ 1 ] = 0.5;
			linearInitializationAnchor[ 2 ] = 0.5;
			linearIterations = 50;
			maxDisplacementRatio = 0.4;
			deformableLevels = 3;
			numberOfFixedImages = 0;
			numberOfRANSACIterations = 5000;
			printLinear = false;
			printStats = false;
			RANSACInlierDistance = 50;
			RANSACMaxScale = 10;
			statIntervalUpdate = 10;
			useRANSAC = true;
			useScale = true;
			writePairs = false;
			writeSingleFileTransforms = false;
			transformSubdirectory = std::string( "transforms" );
			errorMapsSubdirectory = std::string( "errorMaps" );
		};

protected:

	class Measure {

	public:

		float E, landmarkAv, landmarkMax, landmarkSTD;
		Measure( float e ) : landmarkAv( 0 ), landmarkMax( 0 ), landmarkSTD( 0 ),
			E ( e ) {};

	};


	struct Landmark {

		imageIdType image;
		pointIdType point;

	};
	typedef std::vector < Landmark > Landmarks;

	std::vector < Image > images;
	std::vector < Measure > measures;

	void setupStats();

	void setupDeformableTransforms( int level );
	void setupLinearTransforms();

	typedef std::pair< int, vtkMatrix4x4* > RANSACResult;
	void RANSAC( int image );
	RANSACResult RANSACBatch( int image, int nIterations, int batch );

	void transformPoints( const bool &apply = false );

	void updateStats(); // compute Maxwell distribution parameters
	void countInliers(); // displays number of inliers/outliers

	double updateDeformableTransforms( const float alpha ); // returns error value or -1 if diffeomorphism is not guaranteed
	void saveErrorMaps();
	double updateLinearTransforms();

	void displayStats();
	void displayLinearTransforms();

	void saveIndividualDistanceHistograms();
	void saveDistanceHistograms( const char *file );
	void saveMeasures( const char *file );
	void getBoundingBox( vtkBoundingBox &box, const bool &all = false );

	void saveTransforms();
	void saveBoundingBox();

	// reference landmarks
	std::map < std::string, Landmarks > landmarks;
	void computeLandmarkDistances( float e );
	bool saveTransformedLandmarks();
	void saveLandmarkDistances();


	void readAndApplyFixedImagesTransforms();

	picojson::object stats; // to log stats
};
