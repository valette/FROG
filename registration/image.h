#include <vector>

#include <vtkBoundingBox.h>
#include <vtkGeneralTransform.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#include "point.h"

class Image {
public:

	double refTranslation[ 3 ];
	std::vector< Point > points;
	Stats stats;

	void addPoints( vtkBoundingBox &box );

	void transformPoints(  bool apply = false );

	vtkAbstractTransform *transform;
	vtkSmartPointer<vtkImageData> gradient;

	vtkSmartPointer<vtkGeneralTransform> allTransforms;

	Image () : transform( 0 ) {};

};
