#include "image.h"

void Image::transformPoints() {

	for ( auto &point : this->points )
		this->transform->TransformPoint ( point.xyz, point.xyz2 );

}

void Image::expandBoundingBox( float *box ) {

	for ( const auto &point : this->points ) {

		const float *xyz = point.xyz;

		for ( int k = 0; k < 3; k++ ) {

			if ( xyz[ k ] < box[ 2 * k ] ) box[ 2 * k ] = xyz[ k ];
			if ( xyz[ k ] > box[ 2 * k + 1 ] ) box[ 2 * k + 1 ] = xyz[ k ];

		}

	}

}
