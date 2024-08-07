#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>

#include "stats.h"

#define print( v, size ) { for (int __i = 0; __i < size; __i++ ) { std::cout << v[ __i ]; if ( __i < ( size - 1 ) ) std::cout<< " " ;} std::cout << std::endl;}

int Stats::maxSize = 10000;
int Stats::maxIterations = 10000;
float Stats::epsilon = 1e-6;

void Stats::estimateDistribution() {

	const float esp = 1.59576912160573;
	int iteration = 0;
	const float epsilon = Stats::epsilon;

	while ( iteration++ < Stats::maxIterations ) {

		float sum1 = 0;
		float sum2 = 0;
		float sum3 = 0;
		float sum4 = 0;
		float sum5 = 0;

		for ( int i = 0; i < this->size; i++ ) {

			float f1 = ratio * chipdf( samples[ i ] / c1 ) / c1;
			float f2 = ( 1.0 - ratio ) * chipdf( samples[ i ] / c2 ) / c2;
			float t = f1/ ( f1 + f2 + 1e-16);
			float p = this->samples[ i ] * this->weights[ i ];
			sum1 += t * p;
			sum2 += t * this->weights[ i ];
			sum3 += ( 1.0 - t ) * p;
			sum4 += ( 1.0 - t ) * this->weights[ i ];
			sum5 += this->weights[ i ];

		}

		sum2 = std::max( sum2, epsilon);
		sum3 = std::max( sum3, epsilon);
		sum5 = std::max( sum5, epsilon);
		float nc1 = std::max( epsilon, sum1 / sum2 / esp );
		float nc2 = std::max( epsilon, sum3 / sum4 / esp );
		float nRatio = std::max( epsilon, sum2 / sum5 );

		if ( abs( ( c1 - nc1 ) / nc1 ) < 0.001
			&& abs( ( c2 - nc2 ) / nc2 ) < 0.001
			&& abs( ( nRatio - ratio ) / nRatio ) < 0.001 ) {

			c1 = nc1;
			c2 = nc2;
			ratio = nRatio;
			break;

		} else {

			c1 = nc1;
			c2 = nc2;
			ratio = nRatio;

		}

	}

	return;

}

void Stats::displayParameters() {

	int s = this->size;
	std::cout << "c1=" << c1 << ",";
	std::cout << "c2=" << c2 << ",";
	std::cout << "r=" << ratio << ",";
	std::cout << "nSamples=" << s;

	double sum = std::accumulate( samples.begin(), samples.begin() + s , 0.0 );
	double mean = sum / s;

	double sq_sum = std::inner_product( samples.begin(), samples.begin() + s,
		samples.begin(), 0.0 );

	auto max = max_element( samples.begin(), samples.begin() + s );

	double stdev = sqrt( sq_sum / s - mean * mean);
	std::cout <<",max=" << *max
		<< ",mean="<< mean << ",stdev=" << stdev << std::endl;


}

void Stats::saveHistogram( const char *file ) {

	this->getHistogram();

	std::fstream fs;
	fs.open ( file, std::fstream::out | std::fstream::trunc );

	for ( int i = 0; i < this->histogram.size(); i++ )
		fs << histogram[ i ] << std::endl;

	fs.close();

}

void Stats::saveSamples( const char *file, int subsampling ) {

	std::fstream fs;
	fs.open ( file, std::fstream::out | std::fstream::trunc );

	for ( int i = 0; i < this->size; i++ )
		if ( !( i % subsampling ) ) fs << this->samples[ i ] << std::endl;

	fs.close();

}

void Stats::getHistogram( float binSize ) {

	float max = *std::max_element( samples.begin(), samples.begin() + this->size );
	int size = round( max / binSize ) + 1;
	this->histogram.resize( size );
	std::fill( this->histogram.begin(), this->histogram.begin() + size, 0 );

	for ( int i = 0; i < this->size; i++ )
		this->histogram[ round( this->samples[ i ] / binSize ) ]++;

}
