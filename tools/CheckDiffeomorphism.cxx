/*=========================================================================

Program:   CheckDiffeomorphism : Check that a transform is diffeomorphic
Module:    FROG
Language:  C++
Date:      2019/04
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME VolumeTransform
// .SECTION Description

#include <vtkGeneralTransform.h>
#include <vtkImageData.h>
#include <vtkImageResize.h>
#include <vtkMath.h>
#include <vtkTimerLog.h>

#include "../vtkOpenSURF3D/vtkRobustImageReader.h"
#include "readTransform.h"

int main( int argc, char *argv[] ) {

	if ( argc < 3 ) {
		std::cout << "Usage : CheckDiffeomorphism image transform [spacing]" << std::endl;
		exit( 1 );
	}

	char *file = argv[ 1 ];
	vtkGeneralTransform *transform = readTransform( argv[ 2 ] );

	vtkTimerLog *Timer = vtkTimerLog::New();
	vtkRobustImageReader *imageReader = vtkRobustImageReader::New();

	// Load Volume
	std::cout << "load : " << file << std::endl;
	Timer->StartTimer();
	imageReader->SetFileName( file );
	imageReader->Update();
	Timer->StopTimer();
	std::cout << "Image loaded in " << Timer->GetElapsedTime() << "s" << std::endl;
	vtkImageData *image = imageReader->GetOutput();

	if ( argc > 3 ) {

		vtkImageResize *resize = vtkImageResize::New();
		resize->SetResizeMethodToOutputSpacing();
		double sp = atof( argv[ 3 ] );
		cout << "Resizing image with spacing : " << sp << endl;
		resize->SetOutputSpacing( sp, sp, sp );
		resize->SetInputData( image );
		resize->Update();
		image = resize->GetOutput();

	}

	double in[ 3 ], out[ 3 ], derivative[3][3];
	double *bounds = image->GetBounds();
	std::cout << "image bounds :";
	for ( unsigned int i = 0; i < 6; i++) {
		std::cout << bounds[ i ] << " ";
	}
	std::cout << std::endl;

	for ( int i = 0; i < 3; i++ )
		in[ i ] = 0.5* ( bounds[ 2 * i ] + bounds[ 1 + 2 * i ] );

	transform->InternalTransformDerivative( in, out, derivative );
	float d = vtkMath::Determinant3x3( derivative );

	cout << "Jacobian determinant at "
		<< in[ 0 ] << " " << in[ 1 ] << " " << in[ 2 ] << " : " << d << endl;;

	cout << "Derivative : " << endl;
	for ( int i = 0; i < 3; i++ )
		cout << derivative[ i ][ 0 ] << " " << derivative[ i ][ 1 ] << " " << derivative[ i ][ 2 ] << endl;

	double origin[ 3 ], spacing[ 3 ];
	int dimensions[ 3 ];
	vtkIdType inc[ 3 ];
	image->GetOrigin( origin );
	image->GetSpacing( spacing );
	image->GetDimensions( dimensions );
	image->GetIncrements( inc );
	int n = 0;
	cout << "Computing Jacobian determinants..." << endl;

	for ( int k = 0; k < dimensions[ 2 ]; k++ ) {

		int local = 0;

		for ( int j = 0; j < dimensions[ 1 ]; j++ ) {

			for ( int i = 0; i < dimensions[ 0 ]; i++ ) {

				in[ 0 ] = origin[ 0 ] + i * spacing[ 0 ];
				in[ 1 ] = origin[ 1 ] + j * spacing[ 1 ];
				in[ 2 ] = origin[ 2 ] + k * spacing[ 2 ];

				transform->InternalTransformDerivative( in, out, derivative );
				if ( vtkMath::Determinant3x3( derivative ) < 0 ) n++;

			}

		}

	}

	cout << n << " negative jacobian determinant values ("
		<< std::setprecision( 3 ) << (float) 100.0 * n / 
		( dimensions[ 0 ] * dimensions[ 1 ] * dimensions[ 2 ] ) << "%) " << endl;;

}
