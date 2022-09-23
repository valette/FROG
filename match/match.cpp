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

struct CSVrow {
	Descriptor desc;
	vector<float> meta;
};

typedef vector< CSVrow > CSV;

// reads the keypoint list contained in filename (in CSV format)
CSV* readCSVGZ(string filename) {

	std::string line;
	CSV* myCSV = new CSV();
	std::ifstream GZfile( filename, std::ios_base::in | std::ios_base::binary);
	boost::iostreams::filtering_istream file;
	file.push(boost::iostreams::gzip_decompressor());
	file.push(GZfile);

    while(std::getline(file,line)) {

		CSVrow row;
        std::stringstream  lineStream(line);
        std::string        cell;
        while(std::getline(lineStream,cell,',') && (int)cell[0] != 13 ) {

			//cout << cell.length() << " : " << (int)cell[0] << endl;
			if (row.meta.size() < 6) {

				row.meta.push_back(std::stof(cell));

			} else {

				row.desc.push_back(std::stof(cell));

			}

        }

		if (row.meta.size() > 0) myCSV->push_back(row);

    }

	return myCSV;
}

void writeCSV( CSV &csv, const char *fileName) {

	ofstream file;
	file.open(fileName, std::ofstream::out | std::ofstream::trunc);

	for (auto i = 0; i != csv.size(); i++) {

		CSVrow row = csv[ i ];

		for ( auto j = 0; j < row.meta.size(); j++ ) {

			file <<  csv[ i ].meta[ j ] << ",";

		}

		for ( auto j = 0; j < csv[ i ].desc.size(); j++ ) {

			file <<  csv[ i ].desc[ j ];

			if ( j < csv[ i ].desc.size() - 1 ) {

				file << ",";

			} else {

				if ( i < csv.size() ) {

					file << std::endl;

				}

			}

		}

	}

	file.close();

}

// reads the keypoint list contained in filename (in CSV format)
CSV* readCSV(string filename) {

    std::string line;
	CSV* myCSV = new CSV();
	ifstream file( filename );

    while(std::getline(file,line)) {

		CSVrow row;
		std::stringstream  lineStream(line);
		std::string        cell;
		while(std::getline(lineStream,cell,',') && (int)cell[0] != 13 ) {

			if (row.meta.size() < 6) {

				row.meta.push_back(std::stof(cell));

			} else {

				row.desc.push_back(std::stof(cell));
			}

        }

		if (row.meta.size() > 0) myCSV->push_back(row);

    }

	return myCSV;

}

// reads the keypoint list contained in filename (in binary format)
CSV* readBinary(string filename) {

	FILE* file=fopen(filename.c_str(),"rb");
	CSV* myCSV = new CSV();

	while(!feof(file)) {

		CSVrow row;
		float valF;
		int unused;
		unused = fread(&valF, sizeof(float), 1, file);
		row.meta.push_back(valF);
		unused = fread(&valF, sizeof(float), 1, file);
		row.meta.push_back(valF);
		unused = fread(&valF, sizeof(float), 1, file);
		row.meta.push_back(valF);
		unused = fread(&valF, sizeof(float), 1, file);
		row.meta.push_back(valF);
		unused = fread(&valF, sizeof(float), 1, file);
		row.meta.push_back(valF);
		unused = fread(&valF, sizeof(float), 1, file);
		row.meta.push_back(valF);
		row.desc.resize(48);
		unused = fread(row.desc.data(), sizeof(float), 48, file);
		myCSV->push_back(row);

	}

	return myCSV;
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

MatchVect* ComputeMatches(CSV &csv2, CSV &csv1, float threshold, float dist2second, bool matchAll, float anatVal, bool sym = false) {

	MatchVect* matches = new MatchVect();
	float d1, d2;
	int match = 0;
	int end1 = csv1.size();

	for (int i = 0; i < end1 ; i++) {

		d1 = d2 = FLT_MAX;
		int end2 = csv2.size();

		for ( int j = 0; j < end2 ; j++) {

			//Laplacian
			if (csv1[i].meta[4] != csv2[j].meta[4]) continue;

			//Scale
			if ((csv1[i].meta[3]/csv2[j].meta[3] > 1.3) || 
					(csv2[j].meta[3]/csv1[i].meta[3] > 1.3) )
				continue;

			float dist = norm(csv1[i].desc, csv2[j].desc, csv1[i].desc.size() );

			//Anatomical test (euclidian norm after transform)
			if (anatVal != 0){
				float x1 = csv1[i].meta[0];
				float y1 = csv1[i].meta[1];
				float z1 = csv1[i].meta[2];
				float x2 = csv2[j].meta[0];
				float y2 = csv2[j].meta[1];
				float z2 = csv2[j].meta[2];
				
				float euclNorm = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));

				if (euclNorm > anatVal){
					continue;
				}
			}


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

bool compareCSVrow (CSVrow i,CSVrow j) { return (i.meta[5]>j.meta[5]); }

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

		argumentsIndex += 2;
	}

	vector< CSV* > csvs;
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

			CSVrow row;
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
	csvs.resize(nb);

	#pragma omp parallel shared(filenames, csvs, cout)
	{
		//Full list	
		#pragma omp for schedule(dynamic)
		for (auto it = 0 ; it < nb ; ++it) {

			string ext = filenames[it].substr(filenames[it].find_last_of(".") + 1);
			//cout << filenames[it] << endl;
			CSV* mycsv;
			if (  ext == "csv")
				mycsv = readCSV(filenames[it]);
			else if ( ext == "bin")
				mycsv = readBinary(filenames[it]);
			else if ( ext == "gz")
				mycsv = readCSVGZ(filenames[it]);
			else {
				cerr << "Bad file format : " <<  ext << endl;
			}

			float zT = rigids.size() ? rigids[it][2] : 0;

			auto pend = remove_if (mycsv->begin(), mycsv->end(),
				[zT, zmin, zmax] (CSVrow &val){
					float z = val.meta[2] + zT;
					return z < zmin || z > zmax;
				});

			mycsv->erase (pend, mycsv->end());

			#pragma omp critical
			cout << "image " << it << " rigid : "
				<< rigids[it][0] << ", " << rigids[it][1] << ", " << rigids[it][2]
				<<  " before : " << mycsv->size() << " points, after : " << mycsv->size() << endl << flush;

			csvs[it] = mycsv;

		}

	}

	end = std::chrono::system_clock::now();
	cout << " : " << std::chrono::duration<float>(end-start).count() << "s" << endl;
	start = end;
	cout << csvs[ 0 ][ 0 ][ 0 ].desc.size() << " values per descriptor" << endl;
	cout << "Sorting and pruning..." << endl;
	nb = csvs.size();

	#pragma omp parallel shared(filenames, csvs)
	{

		#pragma omp for schedule(dynamic)
		for (auto it = 0 ; it < nb ; ++it) {

			auto rit= remove_if (csvs[it]->begin(), csvs[it]->end(), [sp](CSVrow row){
				return row.meta[5] < sp;
			});

			csvs[it]->erase(rit, csvs[it]->end());
	

			if ( csvs[it]->size() > np) {
				partial_sort(csvs[it]->begin(), csvs[it]->begin()+np, csvs[it]->end(), compareCSVrow);
				csvs[it]->resize(np);
			} /*else {
				sort(csvs[it]->begin(), csvs[it]->end(), compareCSVrow);
			}*/
			cout << ". ("  << csvs[it]->size() << ")"<< flush;

			if ( writePoints ) {
				std::stringstream outfilename;
				outfilename << "points" << it << ".csv";
				cout << " writing " << outfilename.str() << endl;
				writeCSV( *csvs[it], outfilename.str().c_str() );
			}

		}

	}

	end = std::chrono::system_clock::now();
	cout << " : " << std::chrono::duration<float>(end-start).count() << "s" << endl;
	start = end;

	int sum = 0;

	vector< pair<int, int> > indices;
	for (int i = 0 ; i < csvs.size()-1 ; i++) {
		if (target >= 0){
			if (i != target ){
				indices.push_back( make_pair(i, target) );
			}
		} else {
			for (int j = i+1 ; j < csvs.size() ; j++) {
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
	#pragma omp parallel shared(csvs, pairs)
	{
		#pragma omp for reduction(+:sum) schedule(dynamic)
		for (int it = 0 ; it < indices.size() ; it++) {
			MatchVect* matches = ComputeMatches(*csvs[ indices[it].first ], *csvs[ indices[it].second ], dist, dist2second, matchAll, anatVal);
            if (symFlag){
                MatchVect* matchesSym = ComputeMatches(*csvs[ indices[it].second ], *csvs[ indices[it].first ], dist, dist2second, matchAll, anatVal, true);
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

	for (auto it = 0 ; it < csvs.size() ; it++) {

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
		pointIdType nbPoints = csvs[it]->size();

		fwrite(&nbPoints, sizeof(pointIdType), 1, file);

		for (auto rowIt = 0 ; rowIt < csvs[it]->size() ; rowIt++) {

			vector<float>* data = &csvs[it]->at(rowIt).meta;
			fwrite(data->data(), sizeof(float), data->size(), file);

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
