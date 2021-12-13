#include <cmath>
#include <vtkImageCast.h>
#include <vtkImageData.h>
#include <vtkNIFTIImageWriter.h>

#include "../vtkOpenSURF3D/vtkRobustImageReader.h"

int main( int argc, char *argv[] ) {

	vtkRobustImageReader *reader = vtkRobustImageReader::New();
	vtkImageData *average = vtkImageData::New();
	vtkImageData *stdev = vtkImageData::New();
	vtkImageCast *cast = vtkImageCast::New();
	cast->SetOutputScalarTypeToFloat();
	float nImages = argc - 1;

	for ( int i = 1; i < argc; i++ ) {

		char *file = argv[ i ];
		std::cout << "load : " << file << std::endl;
		reader->SetFileName( file );
		reader->Update();
		vtkImageData *image = reader->GetOutput();
		int dims[ 3 ];
		image->GetDimensions( dims );
		int nbVoxels = dims[ 0 ] * dims[ 1 ] * dims[ 2 ];

		if ( i == 1 ) {

			average->DeepCopy( image );
			stdev->DeepCopy( image );
			average->AllocateScalars(VTK_FLOAT, 1);
			stdev->AllocateScalars(VTK_FLOAT, 1);
			float *ptrOut = ( float * ) average->GetScalarPointer();
			float *ptrOut2 = ( float * ) average->GetScalarPointer();

			for ( int j = 0; j < nbVoxels; j++ ) {

				ptrOut[ j ] = 0;
				ptrOut2[ j ] = 0;

			}

		}

		cast->SetInputData( image );
		cast->Update();
		float *ptrIn = ( float * ) cast->GetOutput()->GetScalarPointer();
		float *ptrOut = ( float * ) average->GetScalarPointer();
		float *ptrOut2 = ( float * ) stdev->GetScalarPointer();

		for ( int j = 0; j < nbVoxels; j++ ) {

			ptrOut[ j ] += ptrIn[ j ] / nImages;
			ptrOut2[ j ] += ( ptrIn[ j ] * ptrIn[ j ] ) / nImages;

		}

		image->Delete();

	}

	vtkNIFTIImageWriter *writer = vtkNIFTIImageWriter::New();
	writer->SetInputData( average );
	writer->SetFileName( "average.nii.gz" );
	writer->Write();

	float *ptrOut = ( float * ) average->GetScalarPointer();
	float *ptrOut2 = ( float * ) stdev->GetScalarPointer();
	int dims[ 3 ];
	average->GetDimensions( dims );
	int nbVoxels = dims[ 0 ] * dims[ 1 ] * dims[ 2 ];

	for ( int j = 0; j < nbVoxels; j++ ) {

		ptrOut2[ j ] = sqrt( ptrOut2[ j ] - ptrOut[ j ] * ptrOut[ j ] );

	}

	writer->SetInputData( stdev );
	writer->SetFileName( "stdev.nii.gz" );
	writer->Write();

}
