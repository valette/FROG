#ifndef __readTransform__
#define __readTransform__

#include <string>
#include <limits>

#include <vtkGeneralTransform.h>
#include <vtkBSplineTransform.h>
#include <vtkMatrixToLinearTransform.h>
//#include "vtkMyBSplineTransform.h"

#include "../vtkOpenSURF3D/picojson.h"
using namespace picojson;

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
	transform->Concatenate( linear );

//	std::cout << "translation : " << translation[ 0 ] << " " << translation[ 1 ] << " " << translation[ 2 ] << std::endl;
//	std::cout << "scale : " << scale[ 0 ] << " " << scale[ 1 ] << " " << scale[ 2 ] << std::endl;

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

//		std::cout << "Level " << level++ << std::endl;

//		std::cout << "dims : " << dims[0] << " " << dims[ 1 ] << " " << dims[ 2 ] << std::endl;

		double bboxMin[ 3 ], bboxMax[ 3 ], spacing[ 3 ], origin[ 3 ];
		for (int i = 0 ; i < 3 ; i++) {
			file >> bboxMin[ i ] >> bboxMax[ i ];
//			std::cout << "min : " << bboxMin[ i ] << ", max : " << bboxMax[ i ] << std::endl;
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

//		vtkMyBSplineTransform *trans = vtkMyBSplineTransform::New();
		vtkBSplineTransform *trans = vtkBSplineTransform::New();
		trans->SetCoefficientData( coeffs );
		transform->Concatenate( trans );

		if ( file.eof() ) {
			std::cerr << "Error while loading Transform" << std::endl;
			exit( 1 );
		} 
	}
	file.close();
	return transform;
}

vtkGeneralTransform *readJSON ( const char *file ) {

	ifstream transformFile ( file );
    std::string str, str2;
	while( transformFile ){
		std::getline( transformFile, str2 );
		str.append( str2 );
	}
	transformFile.close();

	cout << str << endl;
	picojson::value v;
	std::string err;
	picojson::parse( v, str.c_str(), str.c_str() + strlen(str.c_str()), &err );
	if ( !err.empty() ) {
	  std::cerr << err << std::endl;
	}
	object trans = v.get<object>();	
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

vtkGeneralTransform *readTransform ( const char *file ) {

	char fin[6];
	char filename[1000];

	strcpy( filename, file );

	if (filename != NULL) {
		char *p;
		for (p = filename; *p; ++p) {
			*p = tolower(*p);
		}
	}

	strcpy ( fin, ".json" );

	if ( strstr( filename, fin) != NULL) {

		return readJSON( file );

	} else {
		return readTFM( file );
	}

}

#endif
