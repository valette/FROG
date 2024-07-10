#include "image.h"

void Image::transformPoints( bool apply ) {

	for ( auto &point : this->points ) {

		this->transform->TransformPoint ( point.xyz, point.xyz2 );
		if ( !apply ) continue;
		for ( int k = 0; k < 3; k++ ) point.xyz[ k ] = point.xyz2[ k ];

	}

}

void Image::addPoints( vtkBoundingBox &box ) {

	for ( const auto &point : this->points ) {
		double p[ 3 ] = { point.xyz[ 0 ], point.xyz[ 1 ], point.xyz[ 2 ] };
		box.AddPoint( p );
	}

}
