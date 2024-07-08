#include <chrono>
#include <stdio.h>

#include "imageGroup.h"

using namespace std;

int main( int argc, char *argv[] ) {

	chrono::time_point< chrono::system_clock> start, end;
	start = chrono::system_clock::now();
	ImageGroup group;

	if ( argc < 2 ) {

		cout << "Usage : frog inputPairs.bin [options]"<< endl;
		cout << "Options : " << endl;

		cout << endl << "*Linear registration:" << endl;
		cout << "-dlinear 0/1  : display linear parameters during registration. default : " << group.printLinear << endl;
		cout << "-lanchor x y z: set initialization anchor relative position. Default : " << 
			group.linearInitializationAnchor[ 0 ] << " " << group.linearInitializationAnchor[ 1 ] <<
			" " << group.linearInitializationAnchor[ 2 ] << endl;
		cout << "-la <value>   : set alpha. Default : " << group.linearAlpha << endl;
		cout << "-li number    : number of iterations. Default : " << group.linearIterations << endl;
		cout << "-s 0/1        : use scale. Default : " << group.useScale << endl;

		cout << endl << "*Deformable registration:" << endl;
		cout << "-da value     : set alpha. Default : " << group.deformableAlpha << endl;
		cout << "-di number    : number of iterations for each level. Default : " << group.deformableIterations << endl;
		cout << "-dl number    : number of levels. Default : " << group.deformableLevels << endl;
		cout << "-g spacing    : initial grid spacing. Default : " << group.initialGridSize << endl;
		cout << "-gd 0/1       : guaranteed diffeomorphism. Default : " << group.guaranteeDiffeomorphism << endl;
		cout << "-gm ratio     : maximal displacement ratio to guarantee diffeomorphism. Default : " << group.maxDisplacementRatio << endl;

		cout << endl << "*EM Weighting:" << endl;
		cout << "-dstats 0/1   : display stats during registration. Default : " << group.printStats << endl;
		cout << "-emi number   : max number of iterations for EM weighting. Default : " << Stats::maxIterations << endl;
		cout << "-si number    : interval update for statistics. Default : " << group.statIntervalUpdate << endl;
		cout << "-se number    : stats epsilon. Default : " << Stats::epsilon << endl;
		cout << "-ss number    : stats maximal sample size. Default : " << Stats::maxSize << endl;
		cout << "-t threshold  : inlier probability threshold. Default : " << group.inlierThreshold << endl;

		cout << endl << "*Registration with fixed images:" << endl;
		cout << "-fi number    : number of fixed images. Default : " << group.numberOfFixedImages << endl;
		cout << "-fd path      : fixed images transforms directory." << endl;
		cout << "-r 0/1        : use RANSAC instead of linear registration. Default : " << group.useRANSAC << endl;
		cout << "-ri number    : number of RANSAC iterations. Default : " << group.numberOfRANSACIterations << endl;
		cout << "-rs maxScale  : maximum allowed scale for RANSAC iterations. Default : " << group.RANSACMaxScale << endl;
		cout << "-rid value    : RANSAC inlier distance. Default : " << group.RANSACInlierDistance << endl;

		cout << endl << "*Measure error with reference landmarks:" << endl;
		cout << "-l path       : path containing reference landmarks." << endl;
		cout << "-il 0/1       : invert landmarks x and y coordinates. Default : " << group.invertLandmarksCoordinates << endl;

		cout << endl << "*Other parameters:" << endl;
		cout << "-nt number    : set number of threads. Default : number of cores" << endl;
		cout << "-mf file      : path+name of measure.csv file." << endl;
		cout << "-wp 0/1       : write pairs, distances and probabilities. Default : " << group.writePairs << endl;
		cout << "-j            : outputs a single big JSON file for each transform. Default : " << group.writeSingleFileTransforms << endl;
		cout << "-ts subdir    : subdirectory where transforms will be written. Default : " << group.transformSubdirectory << endl;

		return( 1 );

	}

	int argumentsIndex = 2;
	char *landmarks = 0;

	while (argumentsIndex < argc) {

		char* key = argv[argumentsIndex];
		char *value = argv[argumentsIndex + 1];

		if ( strcmp( key, "-da" ) == 0) {
			group.deformableAlpha = atof( value );
		}

		if ( strcmp( key, "-dlinear" ) == 0) {
			group.printLinear = atoi( value );
		}

		if ( strcmp( key, "-dstats" ) == 0) {
			group.printStats = atoi( value );
		}

		if ( strcmp( key, "-di" ) == 0 ) {
			group.deformableIterations = atoi( value );
		}

		if ( strcmp( key, "-dl" ) == 0 ) {
			group.deformableLevels = atoi( value );
		}

		if ( strcmp( key, "-emi" ) == 0 ) {
			Stats::maxIterations = atoi( value );
		}

		if ( strcmp( key, "-fi" ) == 0 ) {
			group.numberOfFixedImages = atoi( value );
		}

		if ( strcmp( key, "-fd" ) == 0 ) {
			group.fixedTransformsDirectory = value;
		}

		if ( strcmp( key, "-g" ) == 0 ) {
			group.initialGridSize = atof( value );
		}

		if ( strcmp( key, "-gd" ) == 0 ) {
			group.guaranteeDiffeomorphism = atoi( value );
		}

		if ( strcmp( key, "-gm" ) == 0 ) {
			group.maxDisplacementRatio = atof( value );
		}

		if ( strcmp( key, "-il" ) == 0 ) {
			group.invertLandmarksCoordinates = atoi( value );
		}

		if ( strcmp( key, "-lanchor" ) == 0) {
			for ( int i = 0; i < 3; i++ )
				group.linearInitializationAnchor[ i ] =
					atof( argv[argumentsIndex + 1 + i ] );
		}

		if ( strcmp( key, "-la" ) == 0) {
			group.linearAlpha = atof( value );
		}

		if ( strcmp( key, "-li" ) == 0 ) {
			group.linearIterations = atoi( value );
		}

		if ( strcmp( key, "-nt" ) == 0 ) {
			omp_set_num_threads( atoi( value ) );
		}

		if ( strcmp( key, "-r" ) == 0 ) {
			group.useRANSAC = atoi( value );
		}

		if ( strcmp( key, "-ri" ) == 0 ) {
			group.numberOfRANSACIterations = atoi( value );
		}

		if ( strcmp( key, "-rs" ) == 0 ) {
			group.RANSACMaxScale = atof( value );
		}

		if ( strcmp( key, "-rid" ) == 0 ) {
			group.RANSACInlierDistance = atof( value );
		}

		if ( strcmp( key, "-s" ) == 0 ) {
			group.useScale = atoi( value );
		}

		if ( strcmp( key, "-se" ) == 0 ) {
			Stats::epsilon = atof( value );
		}

		if ( strcmp( key, "-si" ) == 0 ) {
			group.statIntervalUpdate = atoi( value );
		}

		if ( strcmp( key, "-ss" ) == 0 ) {
			Stats::maxSize = atoi( value );
		}

		if ( strcmp( key, "-t" ) == 0 ) {
			group.inlierThreshold = atof( value );
		}

		if ( strcmp( key, "-ts" ) == 0 ) {
			group.transformSubdirectory = std::string( value );
		}

		if ( strcmp( key, "-l" ) == 0 ) {
			landmarks = value;
		}

		if ( strcmp( key, "-mf" ) == 0 ) {
			group.outputFileName = value;
		}

		if ( strcmp( key, "-wp" ) == 0 ) {
			group.writePairs = atoi( value );
		}

		if ( strcmp( key, "-j" ) == 0 ) {
			group.writeSingleFileTransforms = true;
			argumentsIndex--;
		}

		argumentsIndex += 2;

	}

	char *inputFile = argv[1];
	cout << "Reading : " << inputFile << endl;
	group.readPairs( inputFile );
	if ( landmarks ) group.readLandmarks( landmarks );
	group.run();
	end = chrono::system_clock::now();
	cout << "Total time : " << chrono::duration<float>(end-start).count() << "s" << endl;
	return 0;

}
