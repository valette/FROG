#include <stdio.h>
#include <vector>

#include "stats.h"

#define print( v, size ) { for (int __i = 0; __i < size; __i++ ) { std::cout << v[ __i ]; if ( __i < ( size - 1 ) ) std::cout<< " " ;} std::cout << std::endl;}

class Point {

public:

	float original_xyz[ 3 ]; // original coordinates from pair file
	float other[ 3 ]; // scale, laplacian sign, detector response

	float xyz[3]; // current coordinates
	float xyz2[3]; // transformed coordinates

	std::vector< int > linkIds; // link Ids

};
