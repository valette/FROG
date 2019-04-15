#include <vector>

#include <vtkGeneralTransform.h>
#include <vtkImageData.h>

#include "point.h"

class Image {
public:

	double refTranslation[ 3 ];
	vector< Point > points;
	Stats stats;

	void expandBoundingBox( float *box );

	void transformPoints();

	vtkAbstractTransform *transform;
	vtkImageData *gradient;

	vtkGeneralTransform *allTransforms;

	Image () : gradient( 0 ) {};

};
