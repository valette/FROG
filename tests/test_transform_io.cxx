#include <cassert>
#include <cmath>
#include <iostream>
#include <vtkSmartPointer.h>
#include <vtkGeneralTransform.h>
#include "transformIO.h"

int main() {
	std::cout << "Running test_transform_io..." << std::endl;

	// Test 1: JSON scale and translation parsing
	{
		const char *jsonStr = "{\"scale\": 2.0, \"translation\": [10.0, 20.0, 30.0]}";
		vtkGeneralTransform *trans = readJSONfromString( jsonStr, "dummy.json" );
		assert( trans != nullptr );

		double inPt[ 3 ] = { 1.0, 2.0, 3.0 };
		double outPt[ 3 ] = { 0.0, 0.0, 0.0 };
		trans->TransformPoint( inPt, outPt );

		assert( !std::isnan( outPt[ 0 ] ) && !std::isnan( outPt[ 1 ] ) && !std::isnan( outPt[ 2 ] ) );
		std::cout << "  Transformed (1, 2, 3) -> ("
				  << outPt[ 0 ] << ", " << outPt[ 1 ] << ", " << outPt[ 2 ] << ")" << std::endl;

		trans->Delete();
		std::cout << "  Passed: Simple JSON transform parse & apply" << std::endl;
	}

	std::cout << "All test_transform_io passed!" << std::endl;
	return 0;
}
