#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <fstream>
#include <iostream> 
#include <cmath>
#include <omp.h>
#include <assert.h>
#include <cfloat>
#include <stdio.h>
#include <chrono>
#include <thread>

# ifdef USE_SSE_FOR_MATCHING
#include <boost/align/aligned_allocator.hpp>
#include "immintrin.h"
#endif

#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "../tools/pointIdType.h"
#include "../tools/transformIO.h"

namespace fs = boost::filesystem;
using namespace std;
# ifdef USE_SSE_FOR_MATCHING
typedef vector<float, boost::alignment::aligned_allocator<float, 32> > Descriptor;
#else
typedef vector<float> Descriptor;
#endif
typedef pair<pointIdType, pointIdType> Match;
typedef vector<Match> MatchVect;
typedef std::tuple<unsigned short, unsigned short, MatchVect* > pairMatches;
typedef std::array<double, 3> v3;

struct Point {
	Descriptor desc;
	float coordinates[3];
	float transformedCoordinates[3];
	float scale;
	float laplacianSign;
	float response;
};

typedef vector< Point > Points;

// reads the keypoint list contained in filename (in CSV format)
Points* readCSVGZ(string filename) {

	std::string line;
	Points* points = new Points();
	std::ifstream GZfile( filename, std::ios_base::in | std::ios_base::binary);
	boost::iostreams::filtering_istream file;
	file.push(boost::iostreams::gzip_decompressor());
	file.push(GZfile);

    while(std::getline(file,line)) {

		Point row;
        std::stringstream  lineStream(line);
        std::string        cell;
		int count = 0;

		while(std::getline(lineStream,cell,',') && (int)cell[0] != 13 ) {

			float value = std::stof(cell);

			switch( count ) {

				case 0 : row.coordinates[ 0 ] = value; break;
				case 1 : row.coordinates[ 1 ] = value; break;
				case 2 : row.coordinates[ 2 ] = value; break;
				case 3 : row.scale = value; break;
				case 4 : row.laplacianSign = value; break;
				case 5 : row.response = value; break;
				default : row.desc.push_back( value );

			}

			count++;

        }

		if ( count > 6 ) points->push_back( row );

    }

	return points;
}

void writeCSV( Points &points, const char *fileName) {

	ofstream file;
	file.open(fileName, std::ofstream::out | std::ofstream::trunc);

	for (auto i = 0; i != points.size(); i++) {

		Point row = points[ i ];

		file <<  points[ i ].coordinates[ 0 ] << ",";
		file <<  points[ i ].coordinates[ 1 ] << ",";
		file <<  points[ i ].coordinates[ 2 ] << ",";
		file <<  points[ i ].scale << ",";
		file <<  points[ i ].laplacianSign << ",";
		file <<  points[ i ].response << ",";

		for ( auto j = 0; j < points[ i ].desc.size(); j++ ) {

			file <<  points[ i ].desc[ j ];

			if ( j < points[ i ].desc.size() - 1 ) {

				file << ",";

			} else {

				if ( i < points.size() ) {

					file << std::endl;

				}

			}

		}

	}

	file.close();

}

// reads the keypoint list contained in filename (in CSV format)
Points* readCSV(string filename) {

    std::string line;
	Points* points = new Points();
	ifstream file( filename );

    while(std::getline(file,line)) {

		Point row;
		std::stringstream  lineStream(line);
		std::string        cell;
		int count = 0;

		while(std::getline(lineStream,cell,',') && (int)cell[0] != 13 ) {

			float value = std::stof(cell);

			switch( count ) {

				case 0 : row.coordinates[ 0 ] = value; break;
				case 1 : row.coordinates[ 1 ] = value; break;
				case 2 : row.coordinates[ 2 ] = value; break;
				case 3 : row.scale = value; break;
				case 4 : row.laplacianSign = value; break;
				case 5 : row.response = value; break;
				default : row.desc.push_back( value );

			}

			count++;

        }

		if ( count > 6 ) points->push_back( row );

    }

	return points;

}

// reads the keypoint list contained in filename (in binary format)
Points* readBinary(string filename) {

	FILE* file=fopen(filename.c_str(),"rb");
	Points* points = new Points();

	while(!feof(file)) {

		Point row;
		float valF;
		int unused;
		unused = fread(&valF, sizeof(float), 1, file);
		row.coordinates[ 0 ] = valF;
		unused = fread(&valF, sizeof(float), 1, file);
		row.coordinates[ 1 ] = valF;
		unused = fread(&valF, sizeof(float), 1, file);
		row.coordinates[ 2 ] = valF;
		unused = fread(&valF, sizeof(float), 1, file);
		row.scale = valF;
		unused = fread(&valF, sizeof(float), 1, file);
		row.laplacianSign = valF;
		unused = fread(&valF, sizeof(float), 1, file);
		row.response = valF;
		row.desc.resize(48);
		unused = fread(row.desc.data(), sizeof(float), 48, file);
		points->push_back(row);

	}

	return points;
}


# ifdef USE_SSE_FOR_MATCHING

static inline float _mm256_reduce_add_ps(__m256 x) {
    /* ( x3+x7, x2+x6, x1+x5, x0+x4 ) */
    const __m128 x128 = _mm_add_ps(_mm256_extractf128_ps(x, 1), _mm256_castps256_ps128(x));
    /* ( -, -, x1+x3+x5+x7, x0+x2+x4+x6 ) */
    const __m128 x64 = _mm_add_ps(x128, _mm_movehl_ps(x128, x128));
    /* ( -, -, -, x0+x1+x2+x3+x4+x5+x6+x7 ) */
    const __m128 x32 = _mm_add_ss(x64, _mm_shuffle_ps(x64, x64, 0x55));
    /* Conversion to float is a no-op on x86-64 */
    return _mm_cvtss_f32(x32);
}

// 75 op
inline float norm(Descriptor& pts1, Descriptor& pts2, int size) {
	float result = .0;
	float* a = pts1.data();
	float* b = pts2.data();

	for (int i = 0 ; i < size; i+=8) {
		const __m256 x_1Vec = _mm256_load_ps(a+i);
		const __m256 y_1Vec = _mm256_load_ps(b+i);
    		__m256 ans = _mm256_sub_ps(x_1Vec, y_1Vec); // Soustraction
		ans =  _mm256_mul_ps(ans, ans);	// ^2
		result += _mm256_reduce_add_ps(ans);
	}

	return result;
}

#else

inline float norm(Descriptor& pts1, Descriptor& pts2, int size) {
	float result = .0;

	for ( int i = 0 ; i < size; i++ ) {
		result += ( pts1[ i ] - pts2[ i ] ) * ( pts1[ i ] - pts2[ i ] );
	}

	return result;
}

#endif

MatchVect* ComputeMatches(Points &points2, Points &points1, float threshold, float dist2second, bool matchAll, float anatVal, bool sym = false) {

	MatchVect* matches = new MatchVect();
	float d1, d2;
	int match = 0;
	int end1 = points1.size();

	for (int i = 0; i < end1 ; i++) {

		d1 = d2 = FLT_MAX;
		int end2 = points2.size();

		for ( int j = 0; j < end2 ; j++) {

			//Laplacian
			if (points1[i].laplacianSign != points2[j].laplacianSign) continue;

			//Scale
			if ((points1[i].scale/points2[j].scale > 1.3) || 
					(points2[j].scale/points1[i].scale > 1.3) )
				continue;

			//Anatomical test (euclidian norm after transform)
			if (anatVal != 0){
				float x1 = points1[i].transformedCoordinates[0];
				float y1 = points1[i].transformedCoordinates[1];
				float z1 = points1[i].transformedCoordinates[2];
				float x2 = points2[j].transformedCoordinates[0];
				float y2 = points2[j].transformedCoordinates[1];
				float z2 = points2[j].transformedCoordinates[2];
				
				float euclNorm = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));

				if (euclNorm > anatVal){
					continue;
				}
			}

			float dist = norm(points1[i].desc, points2[j].desc, points1[i].desc.size() );

			if (matchAll && sqrt(dist) < threshold) {
				if (sym) {
					matches->push_back(make_pair(i, match));
				} else {
					matches->push_back(make_pair(match, i));
				}
			} else {

				if(dist<d1) {

					d2 = d1;
					d1 = dist;
					match = j;

				} else if(dist<d2) {

					d2 = dist;

				}

			}

		}

		if (!matchAll) {
			if ( ( sqrt(d1/d2) < dist2second || (d2 == FLT_MAX) ) &&
					(sqrt(d1) < threshold)) {

			  if (sym) {
					matches->push_back(make_pair(i, match));
				} else {
					matches->push_back(make_pair(match, i));
				}
			}

		}

	}

	return matches;

}

bool compareCSVrow (Point i,Point j) { return (i.response>j.response); }

int main( int argc, char *argv[] ) {

	std::chrono::time_point<std::chrono::system_clock> start, end;
	int N = 1000000; //something big enough
	float sp = 0; //Seuil pour les points d'interets
	int np = 1000000; //nb point d'interet max
	int nt = std::thread::hardware_concurrency();
	if ( argc < 2 ) {
		std::cout<< "Usage : match pointFiles.txt [options] " << std::endl;
		return 1;
	}
	fs::path full_path = fs::system_complete( fs::path( argv[1] ) );
	float dist = 0.22;
	float dist2second = 1;
	float zmin = -1e20;
	float zmax =  1e20;
	bool matchAll = false;
	bool writePoints = false;
	char* outputFileName = 0;
	int argumentsIndex = 2;
	float anatVal = 0.0;
	bool symFlag = false;
	int target = -1;
	char *transformPrefix = 0;

	while (argumentsIndex < argc) {

		char* key = argv[argumentsIndex];
		char *value = argv[argumentsIndex + 1];

		if (strcmp(key, "-n") == 0) {
			N = atoi(value);
		}

		if (strcmp(key, "-sp") == 0) {
			sp = atof(value);
		}

		if (strcmp(key, "-np") == 0) {
			np = atoi(value);
		}

		if (strcmp(key, "-nt") == 0) {
			nt = atoi(value);
		}

		if (strcmp(key, "-d") == 0) {
			dist = atof(value);
		}

		if (strcmp(key, "-d2") == 0) {
			dist2second = atof(value);
		}

		if (strcmp(key, "-zmin") == 0) {
			zmin = atof(value);
		}

		if (strcmp(key, "-zmax") == 0) {
			zmax = atof(value);
		}

		if (strcmp(key, "-o") == 0) {
			outputFileName = value;
		}

		if (strcmp(key, "-all") == 0) {
			matchAll = true;
		}

		if (strcmp(key, "-p") == 0) {
			writePoints = true;
		}
		
		if (strcmp(key, "-anat") == 0){
			anatVal = atof(value);
		}
		if (strcmp(key, "-sym") == 0){
			symFlag = true;
			argumentsIndex-=1;
        }

		if (strcmp(key, "-targ") == 0){
			target = atoi(value);
		}

		if (strcmp(key, "-transformPrefix") == 0){
			transformPrefix = value;
		}

		argumentsIndex += 2;
	}

	vector< Points* > allPoints;
	fs::directory_iterator end_iter;
	int i = 0;
	vector<v3> rigids;
	std::vector<std::string> filenames;

	if (fs::is_directory(full_path)) {

		for ( fs::directory_iterator dir_itr( full_path );
		      dir_itr != end_iter;
		      ++dir_itr ) {

			if ( fs::is_regular_file( dir_itr->status() ) ) {

				i++;
				filenames.push_back(dir_itr->path().native());

			}

		}

	} else if ( is_regular_file( full_path ) ) {

		std::string line;
		ifstream file( full_path.native() );

		while(std::getline(file,line)) {

			Point row;
		    std::stringstream  lineStream(line);
		    std::string        cell;
			std::getline(lineStream,cell,',');

			if ( cell.find( "/" ) == 0 ) {

				filenames.push_back( cell );
				cout << cell << endl;

			} else {

				filenames.push_back(full_path.parent_path().native() + "/" + cell + ".csv");
				cout << full_path.parent_path().native() + cell << endl;

			}

			v3 point;
			point[ 0 ] = point[ 1 ] = point[ 2 ] = 0;

			try {
				std::getline(lineStream,cell,',');
				point[0] = std::stof (cell);
				std::getline(lineStream,cell,',');
				point[1] = std::stof (cell);
				std::getline(lineStream,cell,',');
				point[2] = std::stof (cell);
			} catch ( ... ) {}

			rigids.push_back(point);

		}

	} else {

		cerr << "Bad argument, first arg must be a valid file or a directory" << endl;
		return 1;

	}

	cout << "Found " << filenames.size() << " files, loading : " << fmin(N, filenames.size()) << endl;
	start = std::chrono::system_clock::now();
	if (filenames.size() > N) filenames.resize(N);
	omp_set_num_threads( nt );
	int nb = filenames.size();
	allPoints.resize(nb);

	#pragma omp parallel shared(filenames, allPoints, cout)
	{
		//Full list	
		#pragma omp for schedule(dynamic)
		for (auto it = 0 ; it < nb ; ++it) {

			string ext = filenames[it].substr(filenames[it].find_last_of(".") + 1);
			//cout << filenames[it] << endl;
			vtkGeneralTransform *transform = 0;
			if ( transformPrefix ) {

				std::string transformFile( transformPrefix );
				transformFile += std::to_string( it );
				transformFile += ".json";
				std::cout << "Reading transform " << transformFile << std::endl;
				transform = readTransform( transformFile.c_str() );

			}

			Points* points;
			if (  ext == "csv")
				points = readCSV(filenames[it]);
			else if ( ext == "bin")
				points = readBinary(filenames[it]);
			else if ( ext == "gz")
				points = readCSVGZ(filenames[it]);
			else {
				cerr << "Bad file format : " <<  ext << endl;
			}

			float zT = rigids.size() ? rigids[it][2] : 0;

			auto pend = remove_if (points->begin(), points->end(),
				[zT, zmin, zmax] (Point &val){
					float z = val.coordinates[2] + zT;
					return z < zmin || z > zmax;
				});

			points->erase (pend, points->end());

			for ( auto point : *points ) {

				if ( transform ) {

					transform->TransformPoint( point.coordinates, point.transformedCoordinates );

				} else {
					for ( int i = 0; i < 3; i++ )
						point.transformedCoordinates[ i ] = point.coordinates[ i ];
				}

			}

			#pragma omp critical
			cout << "image " << it << " rigid : "
				<< rigids[it][0] << ", " << rigids[it][1] << ", " << rigids[it][2]
				<<  " before : " << points->size() << " points, after : " << points->size() << endl << flush;

			allPoints[it] = points;

		}

	}

	end = std::chrono::system_clock::now();
	cout << " : " << std::chrono::duration<float>(end-start).count() << "s" << endl;
	start = end;
	cout << allPoints[ 0 ][ 0 ][ 0 ].desc.size() << " values per descriptor" << endl;
	cout << "Sorting and pruning..." << endl;
	nb = allPoints.size();

	#pragma omp parallel shared(filenames, allPoints)
	{

		#pragma omp for schedule(dynamic)
		for (auto it = 0 ; it < nb ; ++it) {

			auto rit= remove_if (allPoints[it]->begin(), allPoints[it]->end(), [sp](Point row){
				return row.response < sp;
			});

			allPoints[it]->erase(rit, allPoints[it]->end());
	

			if ( allPoints[it]->size() > np) {
				partial_sort(allPoints[it]->begin(), allPoints[it]->begin()+np, allPoints[it]->end(), compareCSVrow);
				allPoints[it]->resize(np);
			} /*else {
				sort(allPoints[it]->begin(), allPoints[it]->end(), compareCSVrow);
			}*/
			cout << ". ("  << allPoints[it]->size() << ")"<< flush;

			if ( writePoints ) {
				std::stringstream outfilename;
				outfilename << "points" << it << ".csv";
				cout << " writing " << outfilename.str() << endl;
				writeCSV( *allPoints[it], outfilename.str().c_str() );
			}

		}

	}

	end = std::chrono::system_clock::now();
	cout << " : " << std::chrono::duration<float>(end-start).count() << "s" << endl;
	start = end;

	int sum = 0;

	vector< pair<int, int> > indices;
	for (int i = 0 ; i < allPoints.size()-1 ; i++) {
		if (target >= 0){
			if (i != target ){
				indices.push_back( make_pair(i, target) );
			}
		} else {
			for (int j = i+1 ; j < allPoints.size() ; j++) {
				indices.push_back( make_pair(i, j) );
			}
		}	
	}

	vector < vector< MatchVect* > > pairs;
	pairs.resize( nb );
	for ( auto i = 0; i < nb; i++ ) {
		pairs[ i ].resize(nb);
		for ( auto j = 0; j < nb; j++ ) pairs[ i ][ j ] = 0;
	}

	cout << "Pairing... " << endl;
	#pragma omp parallel shared(allPoints, pairs)
	{
		#pragma omp for reduction(+:sum) schedule(dynamic)
		for (int it = 0 ; it < indices.size() ; it++) {
			MatchVect* matches = ComputeMatches(*allPoints[ indices[it].first ], *allPoints[ indices[it].second ], dist, dist2second, matchAll, anatVal);
            if (symFlag){
                MatchVect* matchesSym = ComputeMatches(*allPoints[ indices[it].second ], *allPoints[ indices[it].first ], dist, dist2second, matchAll, anatVal, true);
                matches->insert(matches->end(), matchesSym->begin(), matchesSym->end());
            }
			#pragma omp critical
			pairs[ indices[ it ].first ][ indices[ it ].second ] = matches;
			sum += matches->size();
			cout << "." << flush;
		}
	}

	end = std::chrono::system_clock::now();
	cout << " : " << std::chrono::duration<float>(end-start).count() << "s" << endl;
	start = end;

	cout << "Nb Match : " << sum << endl;

	std::stringstream outfilename;

	string f = full_path.stem().string();

	if (outputFileName) {

		string file(outputFileName);
		outfilename << file;

	} else {

		outfilename << "out_" << f << "_" << filenames.size() << ".bin";

	}

	FILE * file = fopen(outfilename.str().c_str(),"wb");

	if (file == NULL) {

		cout << "write error : " << outfilename.str() << endl;
		exit(1);

	}

	unsigned short nbAcq = filenames.size();
	fwrite(&nbAcq, sizeof(unsigned short), 1, file);

	for (auto it = 0 ; it < allPoints.size() ; it++) {

		// Write Filename
		size_t found = filenames[it].find_last_of("/\\");
		string currFile = filenames[it].substr(found+1);
		unsigned short sizeString = currFile.size();
		fwrite(&sizeString, sizeof(unsigned short), 1, file);
		fwrite(currFile.c_str(), sizeof(char), currFile.size(), file);

		// Write Rigid
		v3 tmp;
		if (rigids.size() > 0) {

			tmp = rigids[it];

		} else {

			tmp[0] = 0.0;	tmp[1] = 0.0;	tmp[2] = 0.0;

		}
		
		fwrite(tmp.data(), sizeof(double), 3, file);

		// Write Points
		pointIdType nbPoints = allPoints[it]->size();

		fwrite(&nbPoints, sizeof(pointIdType), 1, file);

		for (auto rowIt = 0 ; rowIt < allPoints[it]->size() ; rowIt++) {

			auto point = &allPoints[it]->at(rowIt);
			fwrite(point->coordinates, sizeof(float), 3, file);
			fwrite(&point->scale, sizeof(float), 1, file);
			fwrite(&point->laplacianSign, sizeof(float), 1, file);
			fwrite(&point->response, sizeof(float), 1, file);

		}

	}

	for ( auto i = 0; i < nb; i++ ) {

		for ( auto j = 0; j < nb; j++ ) {

			// Write local header
			MatchVect* matches = pairs[i][j];
			if ( !matches ) continue;
			unsigned int size = matches->size();
			fwrite(&i, sizeof(unsigned short), 1, file);
			fwrite(&j, sizeof(unsigned short), 1, file);
			fwrite(&size, sizeof(unsigned int), 1, file);
			fwrite(matches->data(), sizeof(Match), size, file);

		}

	}

	fclose(file);
	cout << "Output file : " << outfilename.str() << endl;

}
