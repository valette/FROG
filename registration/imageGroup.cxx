#include <algorithm>
#include <assert.h> 
#include <filesystem>
#include <iostream>
#include <numeric>
#include <omp.h>
#include <sstream>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <vtkBSplineTransform.h>
#include <vtkImageData.h>
#include <vtkLandmarkTransform.h>
#include <vtkMath.h>
#include <vtkMatrixToLinearTransform.h>
#include <vtkPoints.h>

#include "../tools/transformIO.h"

#include "imageGroup.h"

using namespace std;

void checkNaN( float value ) {

	if ( !isnan( value ) ) return;
	cout << endl << "Error : NaN" << endl;
	exit( 1 );

}

void ImageGroup::run() {

	// read transform if there are fixed images
	if ( this->numberOfFixedImages ) readAndApplyFixedImagesTransforms();

	this->setupStats();
	this->setupLinearTransforms();
	this->transformPoints();

	if ( this->useRANSAC && this->numberOfFixedImages ) {

		for ( int i = this->numberOfFixedImages; i < this->images.size(); i++ )
			this->RANSAC( i );

		this->transformPoints();
		this->updateStats();
		if ( this->printStats ) displayStats();
		if ( this->printLinear ) this->displayLinearTransforms();

	} else {

		// linear registration
		cout << endl << "Linear registration" << endl;
		for ( int iteration = 0; iteration < this->linearIterations; iteration++ ) {

			Measure measure;
			cout << "Linear registration, iteration " << iteration + 1
				<< "/" << this->linearIterations << endl;

			if ( !( iteration % this->statIntervalUpdate ) ) this->updateStats();
			if ( this->printStats ) displayStats();
			measure.E = this->updateLinearTransforms();
			if ( this->printLinear ) this->displayLinearTransforms();
			this->transformPoints();
			cout << "E = " << measure.E;
			checkNaN( measure.E );
			if ( !this->computeLandmarkDistances( measure ) ) cout << endl;
			this->measures.push_back( measure );

		}

	}

	this->transformPoints( true );
	this->saveDistanceHistograms( "histograms_linear.csv" );

	if ( this->deformableLevels ) {
		cout << endl << "Deformable registration" << endl;
		vector< int > nGrids;
		this->countInliers();

		for ( int level = 0; level < this->deformableLevels; level++ ) {

			cout << endl << "Level " << level + 1 << "/" << this->deformableLevels << endl;
			this->setupDeformableTransforms( level );
			this->transformPoints();
			Measure measure;
			int numberOfGrids = 1;
			float alpha = this->deformableAlpha;
			cout << "alpha = " << alpha << endl;
			int numberOfDiffeomorphicIterations = 0;

			for ( int iteration = 0; iteration < this->deformableIterations; iteration++ ) {

				cout << "Level " << level + 1 << "/" << this->deformableLevels
					<< ", Iteration " << iteration + 1 <<"/"
					<< this->deformableIterations << endl;

				if ( !( iteration % this->statIntervalUpdate ) ) this->updateStats();
				if ( this->printStats ) this->displayStats();
				measure.E = this->updateDeformableTransforms( alpha );
				cout << "E = " << measure.E;
				checkNaN( measure.E );
				if ( measure.E < 0 ) {

					if ( numberOfDiffeomorphicIterations == 0 ) {

						alpha /= 2;
						cout << "Halving alpha. New Value : " << alpha << endl;

					}

					cout << " creating new grid" << endl;
					numberOfGrids++;
					iteration--;
					this->transformPoints( true );
					this->setupDeformableTransforms( level );
					this->transformPoints();
					numberOfDiffeomorphicIterations = 0;
					continue;

				}

				numberOfDiffeomorphicIterations++;
				this->transformPoints();
				if ( !this->computeLandmarkDistances( measure ) ) cout << endl;
				this->measures.push_back( measure );

			}

			this->countInliers();
			cout << "Number of grids for this level : " << numberOfGrids << endl;
			nGrids.push_back( numberOfGrids );
			this->transformPoints( true );

		}

		int totalNumberOfGrids = 0;
		cout << "Grids per level : ";

		for ( int i = 0; i < nGrids.size(); i++ ) {

			totalNumberOfGrids +=nGrids[ i ];
			cout << nGrids[ i ] << " ";

		}

		cout << endl << "Total number of grids : " << totalNumberOfGrids << endl;
	}

	this->displayStats();
	this->saveDistanceHistograms( "histograms.csv" );
	this->saveMeasures( this->outputFileName );
	this->saveTransforms();
	this->saveBoundingBox();
	this->saveLandmarkDistances();
	this->saveTransformedLandmarks();
	if ( this->writePairs ) this->writeLinksDistances();
	std::fstream fs;
	fs.open( "bbox.json", fstream::out | fstream::trunc );
	fs << picojson::value( this->stats ).serialize();
	fs.close();

}

void ImageGroup::setupDeformableTransforms( int level ) {

	double size = this->initialGridSize / pow( 2, level ); // grid size
	float box[ 6 ];
	double spacing[ 3 ];
	double length[ 3 ];
	int dims[ 3 ];
	double origin[ 3 ];
	this->getBoundingBox( box );

	// add margin to bounding box
	for ( int k = 0; k < 3; k++ ) {

		float length = box[ 1 + 2 * k ] - box[ 2 * k ];
		float margin = length * this->boundingBoxMargin;
		box[ 1 + 2 * k ] += margin;
		box[ 2 * k ] -= margin;

	}

	for ( int k = 0; k < 3; k++ ) {

		length[ k ] = box[ 1 + 2 * k ] - box[ 2 * k ];
		dims[ k ] = round( length[ k ] / size );
		if ( dims[ k ] < 1 ) dims[ k ] = 1;
		spacing[ k ] = length[ k ] / dims[ k ];
		origin[ k ] = box[ 2 * k ] - spacing[ k ];
		dims[ k ] += 3;

	}

	cout << "Bounding box : ";
	print( box, 6 );
	cout << "Box length : ";
	print( length, 3 );
	cout << "Grid origin : ";
	print( origin, 3 );
	cout << "Grid spacing : ";
	print( spacing, 3 );
	cout << "Grid dimensions (control points): ";
	print( dims, 3 );

	#pragma omp parallel for
	for ( int i = this->numberOfFixedImages; i < this->images.size(); i++ ) {

		vtkImageData *coeffs = vtkImageData::New();
		coeffs->SetOrigin( origin );
		coeffs->SetSpacing( spacing );
		coeffs->SetDimensions( dims );
		coeffs->AllocateScalars( VTK_FLOAT, 3 );
		float *p = ( float *) coeffs->GetScalarPointer();

		for ( int i = 0; i < 3 * dims[ 0 ] * dims[ 1 ] * dims[ 2 ]; i++ )
			p[ i ] = 0;

		vtkBSplineTransform *transform = vtkBSplineTransform::New();
		transform->SetCoefficientData( coeffs );
		transform->SetBorderModeToZero();
		Image *image = &this->images[ i ];
		image->allTransforms->Concatenate( transform );
		image->transform = transform;
		if ( image->gradient ) image->gradient->Delete();
		vtkImageData *gradient = vtkImageData::New();
		gradient->SetOrigin( origin );
		gradient->SetSpacing( spacing );
		gradient->SetDimensions( dims );
		gradient->AllocateScalars( VTK_FLOAT, 4 );
		image->gradient = gradient;

	}

}

// inspired from vtkBSplineTransform.cxx
inline void vtkBSplineTransformWeights( double F[ 4 ], double f ) {

  const double sixth = 1.0/6.0;
  const double half = 0.5;

  double f2 = f * f;

  F[ 3 ] = f2 * f * sixth;
  F[ 0 ] = ( f2 - f ) * half - F[ 3 ] + sixth;
  F[ 2 ] = f + F[ 0 ] - F[ 3 ] * 2;
  F[ 1 ] = 1 - F[ 0 ] - F[ 2 ] - F[ 3 ];

}

double ImageGroup::updateDeformableTransforms( const float alpha ) {

	double sDistances = 0, sWeights = 0;

	#pragma omp parallel for reduction( +:sDistances, sWeights )
	for ( int image1 = this->numberOfFixedImages; image1 < this->images.size(); image1++ ) {

		Image &image = this->images[ image1 ];
		float *gradient = ( float* ) image.gradient->GetScalarPointer();
		int dims[ 3 ];
		vtkIdType increments[ 3 ];
		double spacing[ 3 ];
		double origin[ 3 ];
		image.gradient->GetSpacing( spacing );
		image.gradient->GetOrigin( origin );
		image.gradient->GetDimensions( dims );
		image.gradient->GetIncrements( increments );
		double weights[ 3 ][ 4 ];
		int nValues = 4 * dims[ 0 ] * dims[ 1 ] * dims[ 2 ];
		for ( int i = 0; i < nValues; i++ ) gradient[ i ] = 0;
		Stats *statsA = &image.stats;

		for ( auto const &point : image.points ) {

			const float *pos = point.xyz;
			const float *pA = point.xyz2;
			float sWeight = 0;
			float sDisp[ 3 ] = { 0, 0, 0 };
			// compute point displacement

			for ( auto const &link : point.links ) {

				Image *image2 = &this->images[ link.image ];
				Point *pointB = &image2->points[ link.point ];
				float *pB = pointB->xyz2;
				float dist = sqrt( vtkMath::Distance2BetweenPoints( pA, pB ) );
				float probA = statsA->getInlierProbability( dist );
				float probB = image2->stats.getInlierProbability( dist );
				float weight = min( probA, probB );

				sDistances += weight *weight * dist * dist;
				sWeights += weight * weight;

				if ( weight < this->inlierThreshold ) continue;
				weight *= weight;

				for ( int k = 0; k < 3; k++ )
					sDisp[ k ] += weight * ( pB[ k ] - pA[ k ] );

				sWeight += weight;

			}

			for ( auto const &link : point.hardLinks ) {

				Image *image2 = &this->images[ link.image ];
				Point *pointB = &image2->points[ link.point ];
				float *pB = pointB->xyz2;
				float dist = sqrt( vtkMath::Distance2BetweenPoints( pA, pB ) );
				float weight = this->landmarksConstraintsWeight;

				sDistances += weight *weight * dist * dist;
				sWeights += weight * weight;
				weight *= weight;

				for ( int k = 0; k < 3; k++ )
					sDisp[ k ] += weight * ( pB[ k ] - pA[ k ] );

				sWeight += weight;

			}

			if ( sWeight == 0 ) continue;
			// add point displacement to gradient
			int idZ = 0;

			for ( int k = 0; k < 3; k++) {

				float coord = ( pos[ k ] - origin[ k ] ) / spacing[ k ];
				int intCoord = floor( coord );
				vtkBSplineTransformWeights( weights[ k ], coord - intCoord );
				idZ += ( intCoord - 1 ) * increments[ k ];

			}

			for ( int k = 0; k < 4; k++ ) {

				int idY = idZ;

				for ( int j = 0; j < 4; j++ ) {

					int idX = idY;

					for ( int i = 0; i < 4; i++ ) {

						double w = weights[ 0 ][ i ]* weights[ 1 ][ j ] * weights[ 2 ][ k ];

						for ( int l = 0; l < 3; l++ )
							gradient[ idX + l ] +=  w * sDisp[ l ];

						gradient[ idX + 3 ] += w * sWeight;
						idX += increments[ 0 ];

					}

					idY += increments[ 1 ];

				}

				idZ += increments[ 2 ];

			}

		}

		// apply gradient to transform
		float *coeffs = ( float *) ( ( vtkBSplineTransform * ) image.transform )
			->GetCoefficientData()->GetScalarPointer();

		int gradientId = 0;
		int coeffsId = 0;

		for ( int i = 0; i < dims[ 0 ] * dims[ 1 ] * dims[ 2 ]; i++ ) {

			const float coeff = gradient[ gradientId + 3 ];

			if ( coeff > 0 ) {

				for ( int j = 0; j < 3; j++ ) {

					gradient[ gradientId + j ] = coeffs[ coeffsId + j ]
						+ alpha * gradient[ gradientId + j ] / coeff;

				}

			} else {

				for ( int j = 0; j < 3; j++ ) {

					gradient[ gradientId + j ] = coeffs[ coeffsId + j ];

				}

			}

			gradientId += 4;
			coeffsId += 3;

		}

	}

	// substract average displacement on all transforms;
	int nImages = this->images.size();
	float *coeffs[ nImages ];
	int dims[ 3 ];
	double spacing[ 3 ];

	for ( int image = this->numberOfFixedImages; image < nImages; image++ ) {

		vtkImageData *gradient = this->images[ image ].gradient;
		gradient->GetDimensions( dims );
		gradient->GetSpacing( spacing );
		coeffs[ image ] = ( float * ) gradient->GetScalarPointer();

	}

	int nValues = dims[ 0 ] * dims[ 1 ] * dims[ 2 ];
	int nBigCoeffs = 0;
	int direction = 0;
	const float maxD = this->maxDisplacementRatio;
	const bool apply = this->numberOfFixedImages == 0;

	#pragma omp parallel for reduction( + : nBigCoeffs )
	for ( int i = 0; i < nValues; i++ ) {

		int offset = i * 4;

		for ( int j = 0; j < 3; j++ ) {

			double sum = 0;

			if ( apply ) {

				for ( int image = this->numberOfFixedImages; image < nImages; image++ ) {

					sum += coeffs[ image ][ offset + j ];

				}

				sum /= nImages;

			}

			for ( int image = this->numberOfFixedImages; image < nImages; image++ ) {

				coeffs[ image ][ offset + j ] -= sum;

				if ( fabs( coeffs[ image ][ offset + j ] ) > maxD * spacing[ j ] )
					nBigCoeffs++;

			}

		}

	}

	if ( this->guaranteeDiffeomorphism && ( nBigCoeffs > 0 ) ) {

		cout << endl<< "Diffeomorphism is not guaranteed : Iteration canceled" << endl;
		return -1;

	}

	#pragma omp parallel for
	for ( int image1 = this->numberOfFixedImages; image1 < nImages; image1++ ) {

		Image &image = this->images[ image1 ];
		vtkBSplineTransform *transform = ( vtkBSplineTransform * ) image.transform;
		vtkImageData *imageData = transform->GetCoefficientData();
		float *disp = ( float * ) imageData->GetScalarPointer();
		float *newDisp = coeffs[ image1 ];

		int offset = 0;
		int offset2 = 0;

		for ( int i = 0; i < nValues; i++ ) {

			for ( int j = 0; j < 3; j++ ) {

				disp[ offset + j ] = newDisp[ offset2 + j ];

			}

			offset += 3;
			offset2 += 4;

		}

		image.transform->Update();

	}

	return sqrt( sDistances / sWeights );

}

void ImageGroup::updateStats() {

	// add samples to statistics
	#pragma omp parallel for
	for ( int image1 = 0; image1 < this->images.size(); image1++) {

		Image *image = &this->images[ image1 ];
		Stats *stats = &image->stats;
		stats->reset();

		for ( auto const &point : image->points ) {

			const float *pA = point.xyz2;
			for ( auto const &link : point.links ) {

				Image *image2 = &this->images[ link.image ];
				Point *pointB = &image2->points[ link.point ];
				float *pB = pointB->xyz2;
				float dist2 = vtkMath::Distance2BetweenPoints( pA, pB );
				stats->addSample( sqrt( dist2 ) );

			}

		}

		stats->estimateDistribution();

	}

}

void ImageGroup::displayLinearTransforms() {

	for ( int i = this->numberOfFixedImages; i < this->images.size(); i++) {

		cout << "Image " << i << ", translation=";
		auto matrix = ( ( vtkMatrixToLinearTransform *) this->images[ i ].transform )->GetInput( );

		for ( int k = 0; k < 3; k++ ) {

			cout << matrix->GetElement( k, 3 );
			if ( k < 2 ) cout  << " ";

		}

		cout << endl << "scale=";

		for ( int k = 0; k < 3; k++ ) {

			cout << matrix->GetElement( k, k );
			if ( k < 2 ) cout  << " ";

		}

		cout << endl;

	}

}

void ImageGroup::RANSAC( int imageId ) {

	chrono::time_point< chrono::system_clock> start, end;
	start = chrono::system_clock::now();
	cout << "RANSAC registration for image " << imageId << ": "<<std::flush;

	Image *image = &this->images[ imageId ];
	int nBatches = omp_get_num_procs();
	int batchIterations = numberOfRANSACIterations / nBatches;
	vector < RANSACResult > results;

	#pragma omp parallel for
	for ( int i = 0; i < nBatches; i++ ) {

		RANSACResult result = this->RANSACBatch( imageId, batchIterations, i );
		#pragma omp critical
		results.push_back( result );

	}

	vtkMatrixToLinearTransform *transform = ( vtkMatrixToLinearTransform *) image->transform;
	auto matrix = transform->GetInput();
	int maxNumberOfInliers = 0;

	for ( auto res = results.begin(); res != results.end(); res++ ) {

		int nInliers = res->first;

		if ( nInliers > maxNumberOfInliers ) {

			matrix->DeepCopy( res->second );
			maxNumberOfInliers = nInliers;

		}

		res->second->Delete();

	}

	float transformed[ 3 ];
	auto &pts = image->points;
	vtkPoints *source = vtkPoints::New();
	vtkPoints *target = vtkPoints::New();
	source->SetDataTypeToFloat();
	target->SetDataTypeToFloat();
	vtkLandmarkTransform *trans = vtkLandmarkTransform::New();
	trans->SetModeToSimilarity();
	trans->SetSourceLandmarks( source );
	trans->SetTargetLandmarks( target );
	float maxDistance2 = pow( this->RANSACInlierDistance, 2 );

	for ( auto const &pointA : pts ) {

		transform->TransformPoint( pointA.xyz, transformed );

		for ( auto const &link : pointA.links ) {

			Image *image2 = &this->images[ link.image ];
			Point *pointB = &image2->points[ link.point ];
			float *pB = pointB->xyz2;
			if ( vtkMath::Distance2BetweenPoints( transformed, pB ) < maxDistance2 ) {
				source->InsertNextPoint( pointA.xyz );
				target->InsertNextPoint( pB );
			}

		}

	}

	source->Modified();
	target->Modified();
	trans->Update();
	matrix->DeepCopy( trans->GetMatrix() );
	trans->Delete();
	source->Delete();
	target->Delete();
	end = chrono::system_clock::now();
	cout << maxNumberOfInliers << " inliers, computed in " << chrono::duration<float>(end - start).count() << "s" << endl;

	if ( this->stats[ "RANSAC" ].is<picojson::null>() )
		this->stats[ "RANSAC" ] = picojson::value( picojson::array() );

	picojson::object stats;
	stats[ "image" ] = picojson::value( ( double ) imageId );
	stats[ "threshold" ] = picojson::value( this->RANSACInlierDistance );
	stats[ "inliers" ] = picojson::value( (double) maxNumberOfInliers );
	this->stats[ "RANSAC" ].get<picojson::array>().push_back( picojson::value( stats ) );

}

ImageGroup::RANSACResult ImageGroup::RANSACBatch( int imageId, int nIterations, int batch ) {

	int numberOfRansacPoints = 4;
	vtkPoints *source = vtkPoints::New();
	vtkPoints *target = vtkPoints::New();
	source->SetDataTypeToFloat();
	target->SetDataTypeToFloat();
	vtkLandmarkTransform *trans = vtkLandmarkTransform::New();
	trans->SetModeToSimilarity();
	trans->SetSourceLandmarks( source );
	trans->SetTargetLandmarks( target );
	int maxNumberOfInliers = 0;
	float maxDistance2 = pow( this->RANSACInlierDistance, 2 );
	Image *image = &this->images[ imageId ];
	vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
	auto matrix2 = trans->GetMatrix();

	auto &pts = image->points;
	int nPoints = pts.size();
	mt19937 rng( batch * 1000 );

	for ( int i = 0; i < nIterations; i++ ) {

		source->Reset();
		target->Reset();

		for ( int j = 0; j < numberOfRansacPoints; j++ ) {

			while ( true ) {

				int pt = rng() % nPoints;
				Point &point = pts[ pt ];
				auto &links = point.links;
				auto size = links.size();
				if ( size == 0 ) continue;
				source->InsertNextPoint( point.xyz );
				int linkId = rng() % size;
				auto &link = links[ linkId ];
				target->InsertNextPoint(
					this->images[ link.image ].points[ link.point ].xyz );
				break;

			}

		}

		source->Modified();
		target->Modified();
		trans->Update();

		int nInliers = 0;
		float transformed[ 3 ];

		for ( auto const &pointA : pts ) {

			trans->TransformPoint( pointA.xyz, transformed );

			for ( auto const &link : pointA.links ) {

				Image *image2 = &this->images[ link.image ];
				Point *pointB = &image2->points[ link.point ];
				float *pB = pointB->xyz2;
				if ( vtkMath::Distance2BetweenPoints( transformed, pB ) < maxDistance2 )
					nInliers++;

			}

		}

		float determinant = fabs( matrix2->Determinant() );
		if ( ( determinant > this->RANSACMaxScale ) || (  determinant < 1.0 / this->RANSACMaxScale ) ) continue;

		if ( maxNumberOfInliers < nInliers ) {

			maxNumberOfInliers = nInliers;
			matrix->DeepCopy( matrix2 );

		}

	}

	source->Delete();
	target->Delete();
	trans->Delete();
	return make_pair( maxNumberOfInliers, matrix );

}

void ImageGroup::setupLinearTransforms() {

	float box[ 6 ];
	float anchors[ this->images.size() ][ 3 ];
	float averageAnchor[ 3 ] = { 0, 0, 0 };
	float *anchorPosition = this->linearInitializationAnchor;

	for ( int i = 0; i < this->images.size(); i++ ) {

		for ( int j = 0; j < 3; j++ ) {

			box[ 2 * j ] = 1e9;
			box[ 2 * j + 1 ] = -1e9;

		}

		this->images[ i ].expandBoundingBox( box );

		for ( int j = 0; j < 3; j++ ) {

			float c = anchorPosition[ j ];
			float anchor = ( 1 - c ) * box[ 2 * j ] + c * box[ 1 + 2 * j ];
			anchors[ i ][ j ] = anchor;
			if ( i < this->images.size() - this->numberOfFixedImages )
				averageAnchor[ j ] += anchor / ( float ) ( this->images.size() - this->numberOfFixedImages );

		}

	}

	// center volumes
	for ( int i = this->numberOfFixedImages; i < this->images.size() ; i++) {

		Image *image = &this->images[ i ];
		image->transform = vtkMatrixToLinearTransform::New();
		image->allTransforms = vtkGeneralTransform::New();
		image->allTransforms->PostMultiply();
		image->allTransforms->Concatenate( image->transform );

		vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
		matrix->Identity();

		for ( int j = 0; j < 3; j++ )
			matrix->SetElement( j, 3, averageAnchor[ j ] - anchors[ i ][ j ] );

		( ( vtkMatrixToLinearTransform *) image->transform )->SetInput( matrix );

	}

}

void ImageGroup::saveDistanceHistograms( const char *file ) {

	std::fstream fs;
	fs.open( file, fstream::out | fstream::trunc );
	int maxSize = 0;

	for ( int i = 0; i < this->images.size(); i++) {

		Stats *stats = &this->images[ i ].stats;
		stats->getHistogram();
		maxSize = max( maxSize, ( int ) stats->histogram.size() );
		fs << "image " << i ;
		if ( i < this->images.size() - 1 ) fs << ",";
			else fs << endl;

	}

	for ( int d = 0; d < maxSize; d++ ) {

		for ( int i = 0; i < this->images.size(); i++) {

			Stats *stats = &this->images[ i ].stats;

			if ( stats->histogram.size() <= d ) fs << 0;
				else fs << stats->histogram[ d ];

			if ( i < this->images.size() - 1 ) fs << ",";
				else fs << endl;

		}

	}

	fs.close();

}

void ImageGroup::saveIndividualDistanceHistograms() {

	for ( int i = 0; i < this->images.size(); i++) {

		ostringstream file;
		file << "histogram" << i << ".csv";
		this->images[ i ].stats.saveHistogram( file.str().c_str() );

	}

}

void ImageGroup::displayStats() {

	for ( int i = 0; i < this->images.size() ; i++) {

		cout << "Stats " <<  i << ":";
		this->images[ i ].stats.displayParameters();

	}

}

void ImageGroup::transformPoints( bool apply ) {

	#pragma omp parallel for
	for ( int i = this->numberOfFixedImages; i < this->images.size() ; i++ ) {

		auto &image = this->images[ i ];
		image.transformPoints();
		if ( !apply ) continue;
		for ( auto &point : image.points )
			for ( int k = 0; k < 3; k++ ) point.xyz[ k ] = point.xyz2[ k ];

	}

}

bool comparePairs( const vector< float > &p1, const vector< float > & p2) {

    return ( p1[ 4 ] < p2[ 4 ] );

}

void ImageGroup::writeLinksDistances() {

	vector < vector< float > > pairs;

	for ( int i1 = 0; i1 < this->images.size(); i1++ ) {

		Image &image = this->images[ i1 ];
		Stats *stats = &image.stats;

		for ( int p1 = 0; p1 != image.points.size(); p1++ ) {

			auto pointA = &image.points[ p1 ];
			float *pA = pointA->xyz2;

			for ( auto const &link : pointA->links ) {

				int i2 = link.image;
				Image *image2 = &this->images[ i2 ];
				int p2 = link.point;
				Point *pointB = &image2->points[ p2 ];
				float *pB = pointB->xyz2;
				float dist = sqrt( vtkMath::Distance2BetweenPoints( pA, pB ) );
				vector< float > pair;
				pair.push_back( i1 );
				pair.push_back( p1 );
				pair.push_back( i2 );
				pair.push_back( p2 );
				pair.push_back( dist );
				pair.push_back( stats->getInlierProbability( dist ) );
				pairs.push_back( pair );

			}

		}

	}

    sort( pairs.begin(), pairs.end(), comparePairs );

	ofstream file( "pairs.csv.gz", ios_base::out | ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::output> outbuf;
    outbuf.push(boost::iostreams::gzip_compressor());
    outbuf.push(file);
    ostream out(&outbuf);

	for ( int i = 0; i < pairs.size(); i++ ) {

		vector< float > &pair = pairs[ i ];

		for ( int j = 0; j < pair.size(); j++ ) {

			out << pair[ j ];
			if ( j < pair.size() - 1 ) out << ",";

		}

		if ( i  < pairs.size() - 1 ) out << "\n";
	}

	boost::iostreams::close( outbuf );
	file.close();

}

void ImageGroup::countInliers() {

	int nPairs = 0, nInliers = 0, nOutliers = 0;
	picojson::array imagesStats;
	this->stats[ "images" ] = picojson::value( imagesStats );

	// setup json stats record
	for ( int i = 0; i < this->images.size(); i++ ) {

		picojson::object obj;
		this->stats[ "images" ].get<picojson::array>().push_back( picojson::value( obj ) );

	}

	#pragma omp parallel for reduction( +:nPairs, nInliers, nOutliers )
	for ( int image1 = this->numberOfFixedImages; image1 < this->images.size(); image1++ ) {

		Image &image = this->images[ image1 ];
		Stats *statsA = &image.stats;
		int nPairsLocal = 0, nOutliersLocal = 0, nInliersLocal = 0;

		for ( auto const &pointA : image.points ) {

			const float *pA = pointA.xyz2;

			for ( auto const &link : pointA.links ) {

				Image *image2 = &this->images[ link.image ];
				Point *pointB = &image2->points[ link.point ];
				float *pB = pointB->xyz2;
				float dist = sqrt( vtkMath::Distance2BetweenPoints( pA, pB ) );
				float probA = statsA->getInlierProbability( dist );
				float probB = image2->stats.getInlierProbability( dist );
				float weight = min( probA, probB );
				nPairsLocal++;

				if ( weight < this->inlierThreshold )
					nOutliersLocal++;
				else
					nInliersLocal++;

			}

		}

		nPairs += nPairsLocal;
		nInliers += nInliersLocal;
		nOutliers += nOutliersLocal;
		picojson::object stats;
		stats[ "points" ] = picojson::value( ( double ) image.points.size() );
		stats[ "pairs" ] = picojson::value( ( double ) nPairsLocal );
		stats[ "inliers" ] = picojson::value( ( double ) nInliersLocal );
		stats[ "outliers" ] = picojson::value( ( double ) nOutliersLocal );
		picojson::object EMStats;
		EMStats[ "c1" ] = picojson::value( image.stats.c1 );
		EMStats[ "c2" ] = picojson::value( image.stats.c2 );
		EMStats[ "ratio" ] = picojson::value( image.stats.ratio );
		stats[ "EMStats" ] = picojson::value( EMStats );
		this->stats[ "images" ].get<picojson::array>()[ image1 ] = picojson::value( stats );

	}

	cout << "Stats:" << endl;
	cout << nPairs << " half pairs" << endl;
	cout << nInliers << " inliers" << endl;
	cout << nOutliers << " outliers" << endl;
	cout << "Outlier ratio (%): " << ( float ) 100 * nOutliers / nPairs << endl;
	this->stats[ "halfPairs" ] = picojson::value( (double) nPairs );
	this->stats[ "inliers" ] = picojson::value( (double) nInliers );
	this->stats[ "outliers" ] = picojson::value( (double) nOutliers );
	this->stats[ "outlierRatio" ] = picojson::value( (double) nOutliers / nPairs );

}


double ImageGroup::updateLinearTransforms() {

	double sDistances = 0, sWeights = 0;

	#pragma omp parallel for reduction( +:sDistances, sWeights )
	for ( int image1 = this->numberOfFixedImages; image1 < this->images.size(); image1++ ) {

		float diff[3];
		double sDisp[ 3 ] = { 0, 0, 0 };
		double sPosA[ 3 ] = { 0, 0, 0 };
		double sPosA2[ 3 ] = { 0, 0, 0 };
		double sPosB[ 3 ] = { 0, 0, 0 };
		double sPosB2[ 3 ] = { 0, 0, 0 };
		double sWeight = 0;
		Image &image = this->images[ image1 ];
		Stats *statsA = &image.stats;

		for ( auto const &pointA : image.points ) {

			const float *pA = pointA.xyz2;

			for ( auto const &link : pointA.links ) {

				Image *image2 = &this->images[ link.image ];
				Point *pointB = &image2->points[ link.point ];
				float *pB = pointB->xyz2;
				float dist = 0;

				for ( int k = 0; k < 3; k++ ) {

					diff[ k ] = pB[ k ] - pA[ k ];
					dist += diff[ k ] * diff[ k ];

				}

				dist = sqrt( dist );
				float probA = statsA->getInlierProbability( dist );
				float probB = image2->stats.getInlierProbability( dist );
				float weight = min( probA, probB );

				sDistances += weight * weight * dist * dist;
				sWeights += weight * weight;

				for ( int k = 0; k < 3; k++ ) {

					float posA = pA[ k ];
					float posB = pB[ k ];
					sDisp[ k ] += weight * diff[ k ];
					sPosA[ k ] += weight * posA;
					sPosB[ k ] += weight * posB;
					sPosA2[ k ] += weight * posA * posA;
					sPosB2[ k ] += weight * posB * posB;

				}

				sWeight += weight;

			}

		}

		auto matrix = ( ( vtkMatrixToLinearTransform *) image.transform )->GetInput();
		for ( int k = 0; k < 3; k++ ) {

			float sWeight2 = sWeight * sWeight;
			float scale = matrix->GetElement( k, k );
			float newScale = this->useScale ? pow(
				( sWeight * sPosB2[ k ] - sPosB[ k ] * sPosB[ k ] ) /
				( sWeight * sPosA2[ k ] - sPosA[ k ] * sPosA[ k ] ),
				0.5 * this->linearAlpha ) : 1.0;

			if ( isnan( newScale ) ) continue;
			matrix->SetElement( k, k, scale * newScale );
			float translation = matrix->GetElement( k, 3 );
			if ( isnan( translation ) ) continue;
			matrix->SetElement( k, 3, translation
				+ this->linearAlpha * sDisp[ k ] / sWeight
				+ sPosA[ k ] * ( 1 - newScale ) / sWeight );

		}

		image.transform->Update();

	}

	return sqrt( sDistances / sWeights );

}

void ImageGroup::setupStats() {

	// add samples to statistics
	for ( auto &image : this->images )
		for ( auto const &point : image.points )
			for ( auto const &link : point.links )
				image.stats.addSlot();

}

void ImageGroup::addLandmarks( const char *path, bool asConstraints ) {

	vector < string > files;
	std::map < std::string, Landmarks > constraints;

	for( const auto &p: filesystem::directory_iterator( path ) )
		files.push_back ( p.path() );

	sort( files.begin(), files.end() );

	for ( int i = 0; i < files.size(); i++ ) {

		ifstream infile( files[ i ] );
		string line;
		if ( i > this->images.size() - 1 ) continue;

		while ( getline( infile, line ) ) {

			if ( line[ 0 ] == '#' ) continue;
			int pos = line.find( ',');
			string name = line.substr( 0, pos );
			line.erase( 0, pos + 1 );
			Landmark landmark;
			landmark.image = i;
			auto &points = this->images[ i ].points;
			landmark.point = points.size();
			Point pt;

			for ( int j = 0; j < 3; j++ ) {

				int pos = line.find( ',');
				string coord = line.substr( 0, pos );
				line.erase( 0, pos + 1 );				
				pt.xyz[ j ] = stof( coord );
				if ( j < 2 && this->invertLandmarksCoordinates )
					pt.xyz[ j ] *= -1; // get opposite x and y coordinates!
				for ( int k = 0; k < 3; k++ )
					pt.original_xyz[ k ] = pt.xyz[ k ];

			}

			points.push_back( pt );
			this->landmarks[ name ].push_back( landmark );
			constraints[ name ].push_back( landmark );

		}

	}

	if ( !asConstraints ) return;

	for ( auto const &[name, currentLandmarks] : constraints) {

		for ( auto const &landmark : currentLandmarks ) {

			for ( auto const &landmark2 : currentLandmarks ) {

				if ( ( landmark.image == landmark2.image ) &&
					( landmark.point == landmark2.point ) )
					continue;

				this->images[ landmark.image ].points[ landmark.point ].hardLinks.push_back( { landmark2.image, landmark2.point } );

			}

		}
	}

}

bool ImageGroup::computeLandmarkDistances( Measure &measure ) {

	if ( !this->landmarks.size() ) return false;
	vector< float > distances;

	for ( auto const &[name, currentLandmarks] : this->landmarks) {

		float center[ 3 ] = { 0, 0, 0 };

		for ( auto const &landmark : currentLandmarks ) {

			const auto &pt = this->images[ landmark.image ].points[ landmark.point ];
			for ( int k = 0; k < 3; k++ ) center[ k ] += pt.xyz2[ k ] / currentLandmarks.size();

		}

		for ( auto const &landmark : currentLandmarks ) {

			const auto &pt = this->images[ landmark.image ].points[ landmark.point ];
			float d2 = vtkMath::Distance2BetweenPoints( pt.xyz2, center );
			distances.push_back( sqrt( d2 ) );

		}

	}

	double sum = accumulate( distances.begin(), distances.end(), 0.0 );
	double mean = sum / distances.size();

	double sq_sum = std::inner_product( distances.begin(), distances.end(),
		distances.begin(), 0.0 );

	auto max = max_element( distances.begin(), distances.end());

	double stdev = sqrt( sq_sum / distances.size() - mean * mean);
	cout << ", " << distances.size() << " landmarks:max=" << *max
		<< ",average="<< mean << ",stdev=" << stdev << endl;

	measure.landmarkAv = mean;
	measure.landmarkMax = *max;
	measure.landmarkSTD = stdev;
	return true;

}

bool ImageGroup::saveTransformedLandmarks() {

	if ( !this->landmarks.size() ) return false;
	picojson::object transformedLandmarks;

	for ( const auto &[ name, arr ] : this->landmarks ) {

		picojson::array landmarksArray;

		for ( const auto &landmark : arr ) {

			const auto &pt = this->images[ landmark.image ].points[ landmark.point ];

			picojson::object land;
			land[ "image" ] = picojson::value( ( double ) landmark.image );
			picojson::array coords;

			for ( int i = 0; i < 3; i++ )
				coords.push_back( picojson::value( pt.xyz2[ i ] ) );

			land[ "xyz" ] = picojson::value( coords );
			landmarksArray.push_back( picojson::value( land ) );

		}

		transformedLandmarks[ name ] = picojson::value( landmarksArray );

	}

	std::fstream fs;
	fs.open( "transformedLandmarks.json", fstream::out | fstream::trunc );
	fs << picojson::value( transformedLandmarks ).serialize();
	fs.close();
	return true;

}

void ImageGroup::saveLandmarkDistances() {

	if ( !this->landmarks.size() ) return;
	std::fstream fs;
	fs.open( "distances.txt", fstream::out | fstream::trunc );
	picojson::object distances;

	for ( const auto &[ name, landmarks ] : this->landmarks ) {

		float center[ 3 ] = { 0, 0, 0 };

		for ( auto const &landmark : landmarks ) {

			const auto &pt = this->images[ landmark.image ].points[ landmark.point ];
			for ( int k = 0; k < 3; k++ ) center[ k ] += pt.xyz2[ k ] / landmarks.size();

		}

		for ( auto const &landmark : landmarks ) {

			const auto &pt = this->images[ landmark.image ].points[ landmark.point ];
			float d2 = vtkMath::Distance2BetweenPoints( pt.xyz2, center );
			fs << sqrt( d2 ) << "," << name << "," << landmark.image << endl;

		}

	}

	fs.close();

}

void ImageGroup::readPairs( char *inputFile ) {

	FILE *file = fopen( inputFile ,"rb");
	unsigned short nImages;
	int unused;
	unused = fread(&nImages, sizeof(unsigned short), 1, file);
	this->images.resize( nImages );
	int numberOfPairs = 0;

	for ( int i = 0; i < nImages; i++) {

		unsigned short nameLength=18;
		unused = fread( &nameLength, sizeof(unsigned short), 1, file);
		char name[ nameLength + 1 ];
		unused = fread( name, sizeof(char), nameLength, file);
		name[ nameLength ] = 0;
		Image *image = &this->images[ i ];
		unused = fread( image->refTranslation, sizeof( double ), 3, file );
		pointIdType nPoints;
		unused = fread( &nPoints, sizeof( pointIdType ), 1, file );
		image->points.resize( nPoints );

		for ( int j = 0; j < nPoints; j++) {

			Point *pt = & image->points[ j ];
			unused = fread( pt->original_xyz, sizeof( float ), 3, file );
			for ( int k = 0; k < 3; k++ ) pt->xyz[ k ] = pt->original_xyz[ k ];
			unused = fread( pt->other, sizeof( float ), 3, file );

		}

	}

	unsigned short image1, image2;

	while ( unused = fread( &image1, sizeof( unsigned short ), 1, file ) ) {

		unused = fread( &image2, sizeof( unsigned short ), 1, file );
		unsigned int size;
		unused = fread(&size, sizeof(unsigned int), 1, file);

		if ( !size ) {

			cout << "Error : number of pairs is 0" << endl;
			exit( 1 );

		}

		for ( int i = 0; i < size; i++ ) {

			pointIdType p1, p2;
			unused = fread(&p1, sizeof(pointIdType), 1, file);
			unused = fread(&p2, sizeof(pointIdType), 1, file);
			this->images[ image1 ].points[ p1 ].links.push_back( { image2, p2 } );
			this->images[ image2 ].points[ p2 ].links.push_back( { image1, p1 } );

		}

		numberOfPairs += size;

	}

	fclose(file);
	cout << numberOfPairs << " pairs read : " << numberOfPairs * 2 << " half pairs " << endl;

}

void ImageGroup::readAndApplyFixedImagesTransforms() {

	if ( this->fixedTransformsDirectory )
		cout << "Reading transforms in directory " << this->fixedTransformsDirectory << endl;

	#pragma omp parallel for
	for ( int i = 0; i < this->numberOfFixedImages; i++ ) {

		Image *image = &this->images[ i ];

		if ( this->fixedTransformsDirectory ) {

			ostringstream file;
			file << this->fixedTransformsDirectory << "/" << i << ".json";
			image->allTransforms = readTransform( file.str().c_str() );

		} else {

			image->transform = vtkMatrixToLinearTransform::New();
			image->allTransforms = vtkGeneralTransform::New();
			image->allTransforms->PostMultiply();
			image->allTransforms->Concatenate( image->transform );
			vtkMatrix4x4 *matrix = vtkMatrix4x4::New();
			matrix->Identity();
			( ( vtkMatrixToLinearTransform *) image->transform )->SetInput( matrix );

		}

		for ( auto point = image->points.begin(); point != image->points.end(); point++ ) {

			image->allTransforms->TransformPoint( point->xyz, point->xyz2 );

			for ( int k = 0; k < 3; k++ )
				point->xyz[ k ] = point->xyz2[ k ];

		}

	}

}

void ImageGroup::saveTransforms() {

	filesystem::create_directory( this->transformSubdirectory.c_str() );

	// output to .json
	#pragma omp parallel for
	for ( int image1 = this->numberOfFixedImages; image1 < this->images.size(); image1++) {

		vtkGeneralTransform *trans = this->images[ image1 ].allTransforms;
		ostringstream file;
		file << this->transformSubdirectory << "/" << image1 << ".json";
		writeFrogJSON( trans, file.str().c_str(), !this->writeSingleFileTransforms );

	}

}

void ImageGroup::saveMeasures( const char *file ) {

	std::fstream fs;
	fs.open( file, fstream::out | fstream::trunc );
	fs << "Iteration, E, landmarkAv, landmarkMax, landmarkSTD" << endl;

	for ( auto i = this->measures.begin(); i != this->measures.end(); i++ ) {

		fs << i - this->measures.begin() << "," << i->E
			<< "," << i->landmarkAv << "," << i->landmarkMax
			<< ',' << i->landmarkSTD << endl;

	}

	fs.close();

}

void ImageGroup::saveBoundingBox() {

	float box[ 6 ];
	picojson::array min, max;
	this->getBoundingBox( box, true );

	for ( int i = 0; i < 3; i++ ) {

		min.push_back( picojson::value( box[ 2 * i ] ) );
		max.push_back( picojson::value( box[ 1 + 2 * i ] ) );

	}

	picojson::array bbox;
	bbox.push_back( picojson::value( min ) );
	bbox.push_back( picojson::value( max ) );
	this->stats[ "bbox" ] = picojson::value( bbox );

}

void ImageGroup::getBoundingBox( float *box, bool all ) {

	for ( int i = 0; i < 3; i++ ) {

		box[ 2 * i ] = 1e9;
		box[ 2 * i + 1 ] = -1e9;

	}

	for (int i = all ? 0 : this->numberOfFixedImages; i < this->images.size(); i++ ) {

		this->images[ i ].expandBoundingBox( box );

	}

}
