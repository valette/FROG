/*=========================================================================

Program:   MeshTransform : Transform a mesh
Module:    FROG
Language:  C++
Date:      2021/06
Auteur:   Sebastien Valette

=========================================================================*/
// .NAME MeshTransform
// .SECTION Description

#include <vtkGeneralTransform.h>
#include <vtkOBJReader.h>
#include <vtkOBJWriter.h>
#include <vtkPLYReader.h>
#include <vtkTimerLog.h>
#include <vtkTransformPolyDataFilter.h>

#include "transformIO.h"

int main( int argc, char *argv[] ) {
	if ( argc < 3 ) {
		std::cout << "Usage : MeshTransform source [-t transform] [-ti inverse_transform] [-o outputFileName]" << std::endl;
		exit( 1 );
	}

	char *outputFile = 0;

	vtkGeneralTransform *transform = vtkGeneralTransform::New();
	transform->Identity();

	int argumentsIndex = 2;
	vtkTimerLog *timer = vtkTimerLog::New();
	timer->StartTimer();

	while ( argumentsIndex < argc ) {

		char * key = argv[ argumentsIndex ];
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


	vtkTransformPolyDataFilter *polyDataTransform = vtkTransformPolyDataFilter::New();
	polyDataTransform->SetTransform( transform );
	vtkOBJReader *reader = vtkOBJReader::New();
//	vtkPLYReader *reader = vtkPLYReader::New();
	reader->SetFileName( argv[ 1 ] );
	reader->Update();
	
	polyDataTransform->SetInputData( reader->GetOutput() );
	polyDataTransform->Update();
	vtkOBJWriter *writer = vtkOBJWriter::New();
	writer->SetInputData( polyDataTransform->GetOutput() );
	writer->SetFileName( "output.obj" );
	writer->Write();

	// reslice
	timer->StopTimer();
	std::cout << "Transform computed in " << timer->GetElapsedTime() << "s" << std::endl;


}
