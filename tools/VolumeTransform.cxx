/*=========================================================================

Program:   VolumeTransform : Transform a volume
Module:    FROG
Language:  C++
Date:      2013/11
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME VolumeTransform
// .SECTION Description

#include <vtkDataArray.h>
#include <vtkGeneralTransform.h>
#include <vtkImageData.h>
#include <vtkImageReslice.h>
#include <vtkMetaImageWriter.h>
#include <vtkNIFTIImageWriter.h>
#include <vtkPointData.h>
#include <vtkTimerLog.h>
#include <vtkImageFlip.h>

#include "../vtkOpenSURF3D/vtkRobustImageReader.h"
#include "readTransform.h"

int main( int argc, char *argv[] ) {
	if ( argc < 3 ) {
		std::cout << "Usage : VolumeTransform source reference [-t transform] [-ti inverse_transform] [-i interpolation] [-o outputFileName] [-rx reverseX]" << std::endl;
		exit( 1 );
	}

	int reverseX = 0 ;
	char *outputFile = 0;
	int interpolation = 1; // linear

	vtkGeneralTransform *transform = vtkGeneralTransform::New();
	transform->Identity();

	int argumentsIndex = 3;
	while ( argumentsIndex < argc ) {

		char * key = argv[ argumentsIndex ];
		char * value = argv[ argumentsIndex + 1 ];

		if ( strcmp( key ,"-t" ) == 0 ) {
			vtkGeneralTransform *trans2 = readTransform( value );
			trans2->Inverse();
			transform->Concatenate( trans2 );
		}

		if ( strcmp( key ,"-ti" ) == 0 ) {
			transform->Concatenate( readTransform( value ) );
		}

		if ( strcmp( key ,"-o" ) == 0 ) {
			outputFile = value;
		}

		if ( strcmp( key ,"-i" ) == 0 ) {
			interpolation = atoi( value );
		}
		
		if (strcmp(key ,"-rx") == 0) {
			reverseX = atoi(value);
		}		

		argumentsIndex += 2;
	}

	vtkTimerLog *Timer = vtkTimerLog::New();
	vtkImageData *images[ 2 ];
	// If we want to use a reader factory, we have to delete readers ourselves
	vtkRobustImageReader *imageReaders[ 2 ];

	for ( int i = 0; i < 2; i++ ) {

		// Create a reader for image and try to load it
		char *file = argv[ i + 1 ];
		imageReaders[ i ] = vtkRobustImageReader::New();

		// Load Volume
		std::cout << "load : " << file << std::endl;
		Timer->StartTimer();
		imageReaders[ i ]->SetFileName( file );
		imageReaders[ i ]->Update();
		Timer->StopTimer();
		std::cout << "Image loaded in " << Timer->GetElapsedTime() << "s" << std::endl;
		images[ i ] = imageReaders[ i ]->GetOutput();

	}

	double *bounds = images[ 0 ]->GetBounds();
	std::cout << "image bounds :";
	for ( unsigned int i = 0; i < 6; i++) {
		std::cout << bounds[ i ] << " ";
	}
	std::cout << std::endl;

	double center[3], transformedCenter[3];
	for ( unsigned int i = 0; i < 3; i++ ) {
		center[ i ] = 0.5 * ( bounds[ i * 2 ] + bounds[ i * 2 + 1 ] );
	}
	std::cout << "center :" << center[ 0 ] << " " << center[ 1 ] << " " << center[ 2 ] << std::endl;

	transform->TransformPoint( center, transformedCenter );
	std::cout << "transformed center :" << transformedCenter[ 0 ] << " " << transformedCenter[ 1 ] << " " << transformedCenter[ 2 ] << std::endl;

	double valueRange[ 2 ] ;
	images[ 0 ]->GetPointData()->GetScalars()->GetRange( valueRange );

	// reslice
	vtkImageReslice *reslice = vtkImageReslice::New();
	reslice->SetInputData( images[ 0 ] );
	reslice->SetOutputSpacing( images[ 1 ]->GetSpacing() );
	reslice->SetOutputExtent( images[ 1 ]->GetExtent() );
	reslice->SetOutputOrigin( images[ 1 ]->GetOrigin() );
	reslice->SetResliceTransform( transform );
	reslice->SetBackgroundLevel( valueRange[ 0 ] );

	if (interpolation == 0 ) {
		reslice->SetInterpolationModeToNearestNeighbor();
	} else {
		reslice->SetInterpolationModeToLinear();
	}

	Timer->StartTimer();
	vtkObject::GlobalWarningDisplayOff();
	reslice->Update();
	vtkObject::GlobalWarningDisplayOn();
	Timer->StopTimer();
	std::cout << "Transform computed in " << Timer->GetElapsedTime() << "s" << std::endl;

	Timer->StartTimer();

	// write output
	vtkImageWriter *writer;

	if ( outputFile ) {

		char fin[ 10 ];
		char filename[ 1000 ];

		strcpy( filename, outputFile );

		if (filename != NULL) {
			char *p;
			for (p = filename; *p; ++p) {
				*p = tolower(*p);
			}
		}

		strcpy ( fin, ".mhd" );

		if ( strstr( filename, fin) != NULL) {

				writer = vtkMetaImageWriter::New();

		} else {

			strcpy ( fin, ".nii.gz" );
			if ( strstr( filename, fin) != NULL) {

				writer = vtkNIFTIImageWriter::New();

			} else {
				std::cout << "not able to read  " << outputFile << std::endl;
			}

		}

		writer->SetFileName( outputFile );

	} else {

		writer = vtkMetaImageWriter::New();
		writer->SetFileName( "output.mhd" );

	}
	
	vtkImageFlip* flipper = vtkImageFlip::New() ;
    
	if (reverseX) {
        flipper->SetInputConnection(reslice->GetOutputPort()) ;
        flipper->SetFilteredAxis(0) ;
        flipper->Update() ;
        writer->SetInputData(flipper->GetOutput()) ;        
    }
    else
        writer->SetInputData(reslice->GetOutput()) ;
    
	writer->Write() ;
    
	Timer->StopTimer();
	std::cout << "File written in " << Timer->GetElapsedTime() << "s" << std::endl;

   	Timer->Delete();
    flipper->Delete() ;
    reslice->Delete() ;    
    writer->Delete() ;    
    
	for ( int i = 0 ; i < 2 ; ++i )
		imageReaders[ i ]->Delete() ;
}
