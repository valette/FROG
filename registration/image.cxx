#include "image.h"

void Image::transformPoints() {

	Point *point = &this->points[ 0 ];

	for ( int j = 0; j < this->points.size(); j++ ) {

		this->transform->TransformPoint ( point->xyz, point->xyz2 );

		point++;

	}


}

void Image::expandBoundingBox( float *box ) {

	Point *point = &this->points[ 0 ];
	for ( int i = 0; i < this->points.size(); i++) {

		float *xyz = point->xyz;

		for ( int k = 0; k < 3; k++ ) {

			if ( xyz[ k ] < box[ 2 * k ] ) box[ 2 * k ] = xyz[ k ];
			if ( xyz[ k ] > box[ 2 * k + 1 ] ) box[ 2 * k + 1 ] = xyz[ k ];

		}

		point++;
	}

}
