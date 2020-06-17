#include <chrono>
#include <omp.h>
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
		cout << "-da <value>  : set alpha for deformable registration. Default : " << group.deformableAlpha << endl;
		cout << "-dlinear 0/1 : display linear parameters during registration. default : " << group.printLinear << endl;
		cout << "-dstats 0/1  : display stats during registration. Default : " << group.printStats << endl;
		cout << "-di number   : number of iterations for each deformable level. Default : " << group.deformableIterations << endl;
		cout << "-dl number   : number of deformable levels. Default : " << group.deformableLevels << endl;
		cout << "-emi number  : max number of iterations for EM weighting. Default : " << Stats::maxIterations << endl;
		cout << "-g spacing   : initial grid spacing for deformable. Default : " << group.initialGridSize << endl;
		cout << "-gd 0/1      : guaranteed diffeomorphism. Default : " << group.guaranteeDiffeomorphism << endl;
		cout << "-gm ratio    : maximal displacement ratio to guarantee diffeomorphism. Default : " << group.maxDisplacementRatio << endl;
		cout << "-il 0/1      : invert landmarks x and y coordinates. Default : " << group.invertLandmarksCoordinates << endl;
		cout << "-la <value>  : set alpha for linear registration. Default : " << group.linearAlpha << endl;
		cout << "-li number   : number of iterations for linear registration. Default : " << group.linearIterations << endl;
		cout << "-lv <value>  : set linear averaging method. Default : " << group.linearAverageMethod << endl;
		cout << "-nt <number> : set number of threads. Default : number of cores" << endl;
		cout << "-s 0/1       : use scale for linear registration. Default : " << group.useScale << endl;
		cout << "-se number   : stats epsilon. Default : " << Stats::epsilon << endl;
		cout << "-si number   : interval update for statistics. Default : " << group.statIntervalUpdate << endl;
		cout << "-ss number   : stats maximal sample size. Default : " << Stats::maxSize << endl;
		cout << "-t threshold : inlier probability threshold. Default : " << group.inlierThreshold << endl;
		cout << "-l path      : path containing reference landmarks." << endl;
		cout << "-mf measure file      : path+name of measure.csv file." << endl;
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

		if ( strcmp( key, "-la" ) == 0) {
			group.linearAlpha = atof( value );
		}

		if ( strcmp( key, "-li" ) == 0 ) {
			group.linearIterations = atoi( value );
		}

		if ( strcmp( key, "-lv" ) == 0 ) {
			group.linearAverageMethod = atoi( value );
		}

		if ( strcmp( key, "-nt" ) == 0 ) {
			omp_set_num_threads( atoi( value ) );
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

		if ( strcmp( key, "-l" ) == 0 ) {
			landmarks = value;
		}

		if ( strcmp( key, "-mf" ) == 0 ) {
			group.outputFileName = value;
		}

		argumentsIndex += 2;

	}

	if ( landmarks ) group.readLandmarks( landmarks );
	char *inputFile = argv[1];
	cout << "Reading : " << inputFile << endl;
	group.readPairs( inputFile );
	group.run();
	end = chrono::system_clock::now();
	cout << "Total time : " << chrono::duration<float>(end-start).count() << "s" << endl;
	return 0;

}
