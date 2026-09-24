/*=========================================================================

Program:   PointsTransform : Transform a point
Module:    FROG
Language:  C++
Date:      2021/05
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME PointsTransform
// .SECTION Description

#include <fstream>
#include <vtkDataArray.h>
#include <vtkGeneralTransform.h>
#include "transformIO.h"

int main( int argc, char *argv[] ) {

	if ( argc < 3 ) {
		std::cout << "Usage : PointsTransform [-p x y z] [-t transform] [-ti inverse_transform] [-o outputFileName]" << std::endl;
		exit( 1 );
	}

	char *outputFile = 0;
	vtkGeneralTransform *transform = vtkGeneralTransform::New();
	transform->Identity();

	int argumentsIndex = 1;
	double *point = 0;

	while ( argumentsIndex < argc ) {

		char * key = argv[ argumentsIndex ];

		if ( strcmp( key ,"-p" ) == 0 ) {

			if ( argumentsIndex + 3 >= argc ) {
				std::cerr << "Option -p requires 3 coordinates (x y z)." << std::endl;
				return 1;
			}

			point = new double[ 3 ];

			for ( int i = 0; i < 3; i++ ) {

				point[ i ] = atof( argv[ argumentsIndex + 1 + i ] );

			}

			argumentsIndex += 4;
			continue;

		}

		if ( argumentsIndex + 1 >= argc ) {
			std::cerr << "Option " << key << " requires a value." << std::endl;
			return 1;
		}

		char * value = argv[ argumentsIndex + 1 ];

		if ( strcmp( key ,"-t" ) == 0 ) {
			transform->Concatenate( readTransform( value ) );
		}

		if ( strcmp( key ,"-ti" ) == 0 ) {
			vtkGeneralTransform *trans2 = readTransform( value );
			trans2->Inverse();
			transform->Concatenate( trans2 );
		}

		if ( strcmp( key ,"-o" ) == 0 ) {
			outputFile = value;
		}

		argumentsIndex += 2;

	}

	if ( point ) {

		double newPoint[ 3 ];
		cout << "Input point : " << point[ 0 ] << " " << point[ 1 ]
			<< " " << point[ 2 ] << endl;
		transform->TransformPoint( point, newPoint );
		cout << "Output point : " << newPoint[ 0 ] << " " << newPoint[ 1 ]
			<< " " << newPoint[ 2 ] << endl;

		if ( outputFile ) {
			std::ofstream ofs( outputFile );
			if ( ofs.is_open() ) {
				ofs << newPoint[ 0 ] << " " << newPoint[ 1 ] << " " << newPoint[ 2 ] << std::endl;
				ofs.close();
			} else {
				std::cerr << "Could not open output file: " << outputFile << std::endl;
			}
		}


	}


	if ( point ) delete [] point;
	transform->Delete();

}
