#include <math.h>
#include <random>
#include <stdio.h>
#include <vector>

#ifndef STATS_H
#define STATS_H

using namespace std;

// compute Chi probability (from normalized values)
inline float chipdf( float x ) {

	float c =  0.797884560802865;
	float x2 = x*x;
	return c * x2 *exp( - 0.5 * x2 );

}

class Stats {

public:

	static int maxSize;
	static int maxIterations;
	static float epsilon;

	vector< float > samples;
	vector< float > weights;
	int size;
	int virtualSize;
	float c1,c2,ratio;
	mt19937 rng;
	bool needsRandom;

	vector < float > histogram;

	void addSlot() {

		this->virtualSize++;

		if ( this->virtualSize > Stats::maxSize ) {

			this->needsRandom = true;
			return;

		}

		this->samples.push_back( -1 );
		this->weights.push_back( -1 );

	};

	void reset() {

		this->size = 0;

	};

	void addSample( float sample, float weight = 1 ) {

		if ( !this->needsRandom ) {

			this->samples[ this->size ] = sample;
			this->weights[ this->size ] = weight;
			this->size++;
			return;

		}

		if ( this->size == this->samples.size() ) return;
		float random = ( float ) this->rng() / this->rng.max();
		if ( random > (float) this->samples.size() / this->virtualSize  ) return;
		this->samples[ this->size ] = sample;
		this->weights[ this->size ] = weight;
		this->size++;

	};

	void getHistogram( float binSize = 1 );
	void saveHistogram( const char *file );
	void saveSamples( const char *file, int subsampling = 1 );
	void estimateDistribution();
	void displayParameters();

	inline float getInlierProbability( float d ) {

		const float eps = 1e-10;
		float x1 = ratio * chipdf( d / ( c1 + eps ) ) / ( c1 + eps );
		float x2 = ( 1.0 - ratio ) * chipdf( d / ( c2 + eps ) )/ ( c2 + eps);
		return x1 / ( x1 + x2 + eps );

	};

	Stats() : c1( 10 ), c2( 300 ), ratio ( 0.5 ), needsRandom( false ),
		virtualSize( 0 ) {

		this->rng.seed( 0 );

	};

};

#endif
