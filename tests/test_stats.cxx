#include <cassert>
#include <cmath>
#include <iostream>
#include "stats.h"

int main() {
	std::cout << "Running test_stats..." << std::endl;

	// Test 1: Empty stats does not crash or divide by zero
	{
		Stats stats;
		stats.estimateDistribution();
		std::cout << "  Passed: Empty stats estimateDistribution()" << std::endl;
	}

	// Test 2: Inlier probability edge cases
	{
		Stats stats;
		float p0 = stats.getInlierProbability( 0.0f );
		assert( p0 == 1.0f );
		float pSmall = stats.getInlierProbability( 0.05f );
		assert( pSmall == 1.0f );
		float pNormal = stats.getInlierProbability( 10.0f );
		assert( pNormal >= 0.0f && pNormal <= 1.0f );
		std::cout << "  Passed: Inlier probability bounds and edge cases" << std::endl;
	}

	// Test 3: Reservoir sampling
	{
		Stats stats;
		const int numSlots = 10;
		for ( int i = 0; i < numSlots; ++i ) {
			stats.addSlot();
		}
		assert( stats.samples.size() == numSlots );
		assert( stats.size == 0 );

		// Add more samples than slots
		for ( int i = 0; i < 100; ++i ) {
			stats.addSample( static_cast<float>( i ), 1.0f );
		}
		assert( stats.size == numSlots );
		assert( stats.totalSeen == 100 );
		for ( float s : stats.samples ) {
			assert( s >= 0.0f && s < 100.0f );
		}
		std::cout << "  Passed: Reservoir sampling" << std::endl;
	}

	// Test 4: Estimate distribution with valid samples
	{
		Stats stats;
		for ( int i = 0; i < 50; ++i ) {
			stats.addSlot();
		}
		for ( int i = 0; i < 50; ++i ) {
			stats.addSample( 5.0f + ( i % 5 ), 1.0f );
		}
		stats.estimateDistribution();
		assert( !std::isnan( stats.c1 ) && !std::isinf( stats.c1 ) );
		assert( !std::isnan( stats.c2 ) && !std::isinf( stats.c2 ) );
		assert( !std::isnan( stats.ratio ) && !std::isinf( stats.ratio ) );
		std::cout << "  Passed: Distribution estimation" << std::endl;
	}

	std::cout << "All test_stats passed!" << std::endl;
	return 0;
}
