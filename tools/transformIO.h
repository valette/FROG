#ifndef __readTransform__
#define __readTransform__

#include <string>
#include <limits>

#include <vtkGeneralTransform.h>
#include <vtkImageData.h>
#include <vtkBSplineTransform.h>
#include <vtkMatrixToLinearTransform.h>

#include "../vtkOpenSURF3D/picojson.h"
using namespace picojson;

void writeTFM( vtkGeneralTransform *generalTransform, const char *fileName ) {

	fstream fs;
	fs.open( fileName, fstream::out | fstream::trunc );

	vtkAbstractTransform *trans = generalTransform->GetConcatenatedTransform( 0 );
	vtkMatrixToLinearTransform *linear = ( vtkMatrixToLinearTransform *) trans;
	vtkMatrix4x4 *matrix = linear->GetInput();

	for ( int i = 0; i < 3; i++ ) fs << matrix->GetElement( i, 3 ) << " ";
	fs << "-123456 ";

	for ( int i = 0; i < 3; i++ )
		fs << matrix->GetElement( i, i ) << ( i < 2 ? " " : "" );

	fs << endl;

	for ( int i = 1; i < generalTransform->GetNumberOfConcatenatedTransforms(); i++ ) {

		vtkBSplineTransform *transform = ( vtkBSplineTransform * )
			generalTransform->GetConcatenatedTransform( i );

		vtkImageData *imageData = transform->GetCoefficientData();

		int dims[ 3 ];
		double origin[ 3 ];
		double spacing[ 3 ];

		imageData->GetDimensions( dims );
		imageData->GetSpacing( spacing );
		imageData->GetOrigin( origin );

		for ( int k = 0; k < 3; k ++ ) fs << dims[ k ] - 3 << " ";

		for ( int k = 0; k < 3; k++ ) {

			fs << origin[ k ] + spacing[ k ] << " ";
			fs << origin[ k ] + spacing[ k ] * ( dims[ k ] - 2 );
			if ( k < 2 ) fs << " ";
				else fs << endl;

		}

		int count = 0;
		float *values = ( float * ) imageData->GetScalarPointer();
		int nValues = dims[ 0 ] * dims[ 1 ] * dims[ 2 ];

		for ( int j = 0; j < nValues; j++ ) {

			for ( int k = 0; k < 3; k++ ) fs << values[ count++ ] << " ";

			fs << "-123456 -123456" << endl;

		}

	}

	fs.close();

}

void writeFrogJSON( vtkGeneralTransform *generalTransform, const char *fileName ) {

	picojson::array transforms;

	for ( int i = 0; i < generalTransform->GetNumberOfConcatenatedTransforms(); i++ ) {

		auto transform = generalTransform->GetConcatenatedTransform( i );
		auto name = transform->GetClassName();
		picojson::object trans;
		trans[ "type" ] = picojson::value( name );

		if ( strcmp( name, "vtkMatrixToLinearTransform" ) == 0 ) {

			picojson::array matrix;
			auto *m =  ( ( vtkMatrixToLinearTransform *) transform )->GetInput();

			for ( int i = 0; i < 4; i++ ) {

				for ( int j = 0; j < 4; j++ ) {

					matrix.push_back( picojson::value( m->GetElement( i, j ) ) );

				}

			}

			trans[ "matrix" ] = picojson::value( matrix );

		} else if ( strcmp( name, "vtkBSplineTransform" ) == 0 ) {

			vtkImageData *imageData = ( ( vtkBSplineTransform * ) transform )->GetCoefficientData();
			int dims[ 3 ];
			double ori[ 3 ];
			double sp[ 3 ];
			imageData->GetDimensions( dims );
			imageData->GetSpacing( sp );
			imageData->GetOrigin( ori );
			picojson::array dimensions;
			picojson::array origin;
			picojson::array spacing;

			for ( int k = 0; k < 3; k ++ ) {

				dimensions.push_back( picojson::value( ( double ) dims[ k ] ) );
				origin.push_back( picojson::value( ori[ k ] ) );
				spacing.push_back( picojson::value( sp[ k ] ) );

			}

			trans[ "dimensions" ] = picojson::value( dimensions );
			trans[ "origin" ] = picojson::value( origin );
			trans[ "spacing" ] = picojson::value( spacing );
			int count = 0;
			float *values = ( float * ) imageData->GetScalarPointer();
			int nValues = 3 * dims[ 0 ] * dims[ 1 ] * dims[ 2 ];
			picojson::array coeffs;

			for ( int j = 0; j < nValues; j++ ) {

				coeffs.push_back ( picojson::value( values[ j ] ) );

			}

			trans[ "coeffs" ] = picojson::value( coeffs );

		}

		transforms.push_back( picojson::value( trans ) );

	}

	std::fstream fs;
	fs.open( fileName, fstream::out | fstream::trunc );
	picojson::object root;
	root[ "transforms" ] = picojson::value( transforms );
	fs << picojson::value( root ).serialize();
	fs.close();

}

vtkGeneralTransform *readTFM ( const char *fileName ) {

	std::ifstream file( fileName, std::ios::in );

	if( file.fail() ) {

		std::cerr << "Cannot read transform file " << fileName << std::endl;
		exit( 0 );

	}

	double translation[ 3 ];
	file >> translation[ 0 ] >> translation[ 1 ] >> translation[ 2 ];
	double magicNumber;
	file >> magicNumber;
	vtkGeneralTransform *transform = vtkGeneralTransform::New();
	transform->PostMultiply();
	vtkMatrixToLinearTransform *linear = vtkMatrixToLinearTransform::New();
	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	matrix->Identity();

	double scale[ 3 ] = { 1, 1, 1 };

	if ( magicNumber == -123456 ) {

		// we have scales
		file >> scale[ 0 ] >> scale[ 1 ] >> scale[ 2 ];

	}

	for ( int i = 0; i < 3; i++) {

		matrix->SetElement( i, i, scale[ i ] );
		matrix->SetElement( i, 3, translation[ i ] );

	}

	file.ignore ( std::numeric_limits< std::streamsize >::max(), '\n' );
	linear->SetInput( matrix );
	linear->Update();
	transform->Concatenate( linear );

	int level = 0;
	while ( 1 ) {
		if ( file.eof() )
			break;
		bool end = false;
		int dims[ 3 ];
		for (int i = 0; i < 3; i++) {
			float value;
			file >> value;
			if (file.eof()) {
				end = true;
				break;
			}
			dims[ i ] = value;
		}
		if ( end ) break;

		double bboxMin[ 3 ], bboxMax[ 3 ], spacing[ 3 ], origin[ 3 ];

		for (int i = 0 ; i < 3 ; i++) {

			file >> bboxMin[ i ] >> bboxMax[ i ];
			spacing[ i ] = (double) ( bboxMax[ i ] - bboxMin[ i ] ) / dims[ i ];
			dims[ i ] += 3;
			origin[ i ] = bboxMin[ i ] - spacing[ i ];

		}

		file.ignore ( std::numeric_limits<std::streamsize>::max(), '\n' );

		vtkImageData *coeffs = vtkImageData::New();
		coeffs->SetOrigin( origin );
		coeffs->SetSpacing( spacing );
		coeffs->SetDimensions( dims );
		coeffs->AllocateScalars( VTK_FLOAT, 3);
		float *p = ( float *) coeffs->GetScalarPointer();
		int nb = dims[ 0 ] * dims[ 1 ] * dims[ 2 ];

		for ( int i = 0; i < nb ; i++ ) {

			file >> p[ 0 ] >> p[ 1 ] >> p[ 2 ];
			file.ignore ( std::numeric_limits<std::streamsize>::max(), '\n' );	
			p += 3;

		}

		vtkBSplineTransform *trans = vtkBSplineTransform::New();
		trans->SetCoefficientData( coeffs );
		trans->Update();
		transform->Concatenate( trans );

		if ( file.eof() ) {

			std::cerr << "Error while loading Transform" << std::endl;
			exit( 1 );

		} 

	}

	file.close();
	return transform;

}

vtkGeneralTransform *readFrogJSON ( picojson::object root ) {

	vtkGeneralTransform *transform = vtkGeneralTransform::New();
	transform->PostMultiply();

	picojson::array transforms = root[ "transforms" ].get<picojson::array>();

	for ( auto it = transforms.begin(); it != transforms.end(); it++) {

        object trans = it->get<object>();
		std::string type = trans[ "type" ].get<std::string>();

		if ( type.compare ( "vtkMatrixToLinearTransform" ) == 0 ) {

			vtkMatrixToLinearTransform *linear = vtkMatrixToLinearTransform::New();
			vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
			picojson::array m = trans[ "matrix" ].get< picojson::array >();
			int index = 0;

			for ( int i = 0; i < 4; i++ ) {

				for ( int j = 0; j < 4; j++ ) {

					matrix->SetElement( i, j, m[ index++ ].get< double >() );

				}

			}

			linear->SetInput( matrix );
			linear->Update();
			transform->Concatenate( linear );

		} else if ( type.compare ( "vtkBSplineTransform" ) == 0 ) {

			double sp[ 3 ], ori[ 3 ];
			int dims[ 3 ];

			picojson::array dimensions = trans[ "dimensions" ].get< picojson::array >();
			picojson::array origin = trans[ "origin" ].get< picojson::array >();
			picojson::array spacing = trans[ "spacing" ].get< picojson::array >();

			for (int i = 0 ; i < 3 ; i++) {

				dims[ i ] = dimensions[ i ].get< double >();
				sp[ i ] = spacing[ i ].get< double >();
				ori[ i ] = origin[ i ].get< double >();

			}

			vtkImageData *coefficients = vtkImageData::New();
			coefficients->SetOrigin( ori );
			coefficients->SetSpacing( sp );
			coefficients->SetDimensions( dims );
			coefficients->AllocateScalars( VTK_FLOAT, 3);
			float *p = ( float *) coefficients->GetScalarPointer();
			int nb = 3 * dims[ 0 ] * dims[ 1 ] * dims[ 2 ];
			picojson::array coeffs = trans[ "coeffs" ].get< picojson::array >();

			for ( int i = 0; i < nb ; i++ ) {

				p[ i ] = coeffs[ i ].get< double >();

			}

			vtkBSplineTransform *bspline = vtkBSplineTransform::New();
			bspline->SetCoefficientData( coefficients );
			bspline->Update();
			transform->Concatenate( bspline );

		}

	}

	return transform;

}

vtkGeneralTransform *readJSONfromString ( const char *str ) {

	picojson::value v;
	std::string err;
	picojson::parse( v, str, str + strlen( str ) );
	if ( !err.empty() ) std::cerr << err << std::endl;
	object trans = v.get<object>();	

	if ( !( trans[ "transforms" ].is<picojson::null>() ) ) return readFrogJSON( trans );

	cout << str << endl;
	double scale = trans[ "scale" ].get<double>();
	cout << scale <<endl;
	array translation = trans[ "translation" ].get<array>();
	double T[ 3 ];

	for ( int i = 0; i < 3; i++ ) {

		T[ i ] = translation[ i ].get<double>();
		cout << T[ i ] << " ";

	}

	cout << endl;
	vtkGeneralTransform *transform = vtkGeneralTransform::New();
	transform->Scale( scale, scale, scale );
	transform->Translate( T[ 0 ], T[ 1 ], T[ 2 ] );
	return transform;

}

vtkGeneralTransform *readJSON( const char *fileName ) {

	std::ifstream file( fileName );
	std::string str( ( std::istreambuf_iterator< char >( file ) ),
		std::istreambuf_iterator< char >() );

	return readJSONfromString( str.c_str() );

}


vtkGeneralTransform *readTransform( const char *file ) {

	char fin[6];
	char filename[1000];
	strcpy( filename, file );

	if (filename != NULL) {

		char *p;
		for ( p = filename; *p; ++p ) *p = tolower(*p);

	}

	strcpy ( fin, ".json" );
	if ( strstr( filename, fin) != NULL) return readJSON( file );
	else return readTFM( file );

}

#endif
