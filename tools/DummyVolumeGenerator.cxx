#include <vtkImageData.h>
#include <vtkMetaImageWriter.h>
#include <vtkImageWriter.h>

#include "../vtkOpenSURF3D/picojson.h"

using namespace picojson;

int main( int argc, char *argv[] ) {

	if ( argc < 3 ) {

		cout << "Usage : DummyVolumeGenerator bbox.json spacing" << endl;
		exit( -2 );

	}

	char *file = argv[ 1 ];
	double spacing = atof( argv[ 2 ] );

	ifstream bboxFile( file );
    std::string str, str2;

	while( bboxFile ){

		std::getline( bboxFile, str2 );
		str.append( str2 );

	}

	bboxFile.close();
	picojson::value v;
	std::string err;
	picojson::parse( v, str.c_str(), str.c_str() + strlen(str.c_str()), &err );

	if ( !err.empty() ) {

	  std::cerr << err << std::endl;

	}

	object trans = v.get<object>();	
	array bbox = trans[ "bbox" ].get<array>();
	array aMin = bbox[ 0 ].get<array>();
	array aMax = bbox[ 1 ].get<array>();
	double min[ 3 ], max[ 3 ];

	for ( int i = 0; i < 3; i++ ) {

		min[ i ] = aMin[ i ].get<double>();
		max[ i ] = aMax[ i ].get<double>();

	}

	vtkImageData *Volume = vtkImageData::New();
	Volume->SetOrigin( min[ 0 ], min[ 1 ], min[ 2 ] );
	Volume->SetSpacing( spacing, spacing, spacing );
	int dims[ 3 ];

	for ( int i = 0; i < 3; i++ ) {

		dims[ i ] = ceil( ( max[ i ] - min[ i ] ) / spacing );

	}

	Volume->SetDimensions( dims );
	Volume->AllocateScalars(VTK_FLOAT, 1);
	vtkMetaImageWriter *writer =   vtkMetaImageWriter::New();
	writer->SetInputData(Volume);
	writer->SetFileName( "dummy.mhd" );
	writer->Write();
	Volume->Delete();

}
