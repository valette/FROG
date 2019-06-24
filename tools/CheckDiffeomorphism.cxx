/*======================================================================

Program:   CheckDiffeomorphism : Check that a transform is diffeomorphic
Module:    FROG
Language:  C++
Date:      2019/04
Auteur:   Sebastien Valette

======================================================================*/
// .NAME VolumeTransform
// .SECTION Description

#include <vtkGeneralTransform.h>
#include <vtkImageData.h>
#include <vtkImageResize.h>
#include <vtkMath.h>

#include "../vtkOpenSURF3D/vtkRobustImageReader.h"
#include "transformIO.h"

int main( int argc, char *argv[] ) {

	if ( argc < 3 ) {

		std::cout << "Usage : CheckDiffeomorphism image transform [spacing]" << std::endl;
		exit( 1 );

	}

	// Load Volume
	std::cout << "load : " << argv[ 1 ] << std::endl;
	vtkRobustImageReader *imageReader = vtkRobustImageReader::New();
	imageReader->SetFileName( argv[ 1 ] );
	imageReader->Update();
	vtkImageData *image = imageReader->GetOutput();

	vtkGeneralTransform *transform = readTransform( argv[ 2 ] );

	if ( argc > 3 ) {

		double sp = atof( argv[ 3 ] );

		if ( sp > 0 ) {

			cout << "Resizing image with spacing : " << sp << endl;
			vtkImageResize *resize = vtkImageResize::New();
			resize->SetResizeMethodToOutputSpacing();
			resize->SetOutputSpacing( sp, sp, sp );
			resize->SetInputData( image );
			resize->Update();
			image = resize->GetOutput();

		}

	}

	double origin[ 3 ], spacing[ 3 ];
	int dimensions[ 3 ];
	vtkIdType inc[ 3 ];
	image->GetOrigin( origin );
	image->GetSpacing( spacing );
	image->GetDimensions( dimensions );
	image->GetIncrements( inc );
	int n = 0;
	cout << "Computing Jacobian determinants..." << endl;

	#pragma omp parallel for reduction ( + : n )
	for ( int k = 0; k < dimensions[ 2 ]; k++ ) {

		for ( int j = 0; j < dimensions[ 1 ]; j++ ) {

			for ( int i = 0; i < dimensions[ 0 ]; i++ ) {

				double in[ 3 ], out[ 3 ], der[3][3];
				in[ 0 ] = origin[ 0 ] + i * spacing[ 0 ];
				in[ 1 ] = origin[ 1 ] + j * spacing[ 1 ];
				in[ 2 ] = origin[ 2 ] + k * spacing[ 2 ];
				transform->InternalTransformDerivative( in, out, der );
				if ( vtkMath::Determinant3x3( der ) < 0 ) n++;

			}

		}

	}

	cout << n << " negative jacobian determinant values ("
		<< std::setprecision( 3 ) << (float) 100.0 * n / 
		( dimensions[ 0 ] * dimensions[ 1 ] * dimensions[ 2 ] ) << "%) " << endl;;

	return n > 0;

}
