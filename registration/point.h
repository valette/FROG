#include <stdio.h>
#include <vector>

#include "stats.h"
#include "../tools/pointIdType.h"

#define print( v, size ) { for (int __i = 0; __i < size; __i++ ) { std::cout << v[ __i ]; if ( __i < ( size - 1 ) ) std::cout<< " " ;} std::cout << std::endl;}

struct Link {

	unsigned short image;
	pointIdType point;

};


class Point {

public:

	float original_xyz[ 3 ]; // original coordinates from pair file
	float other[ 3 ]; // scale, laplacian sign, detector response

	float xyz[3]; // current coordinates
	float xyz2[3]; // transformed coordinates

	std::vector< Link > links; // links

};
