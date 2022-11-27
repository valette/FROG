import os.path
import vtk
import json

def readTransform( file ):

	transforms = vtk.vtkGeneralTransform()
	transforms.PostMultiply();

	with open( file ) as content:

		data = json.load( content )

		for transform in data[ 'transforms' ]:

			ttype = transform[ 'type' ]

			if ttype == "vtkMatrixToLinearTransform" :

				linear = vtk.vtkMatrixToLinearTransform()
				matrix = vtk.vtkMatrix4x4();
				array = transform[ "matrix" ];
				index = 0;

				for i in range ( 0, 4 ):
					for j in range ( 0, 4 ):
						matrix.SetElement( i, j, array[ index ] );
						index += 1

				linear.SetInput( matrix );
				linear.Update();
				transforms.Concatenate( linear );

			elif ttype == "vtkBSplineTransform" :

				coefficients = vtk.vtkImageData()

				if "file" in transform :

					niiFile = transform[ "file" ]
					reader = vtk.vtkNIFTIImageReader()
					transDir = os.path.dirname( file )
					reader.SetFileName( os.path.join( transDir, niiFile ) )
					reader.Update()
					coefficients.ShallowCopy( reader.GetOutput() )
					origin = [ 0, 0, 0 ]
					qForm = reader.GetQFormMatrix()
					if not qForm : qForm = vtkMatrix4x4()

					for i in range( 3 ) :
						origin[ i ] = qForm.GetElement( i, 3 )

					coefficients.SetOrigin( origin )

				else :

					dimensions = transform[ "dimensions" ]
					origin = transform[ "origin" ]
					spacing = transform[ "spacing" ]

					coefficients.SetOrigin( origin )
					coefficients.SetSpacing( spacing )
					coefficients.SetDimensions( dimensions )
					coefficients.AllocateScalars( vtk.VTK_FLOAT, 3 )
					scalars = coefficients.GetPointData().GetScalars()
					coeffs = transform[ "coeffs" ]
					n = 3 * dimensions[ 0 ] * dimensions[ 1 ] * dimensions[ 2 ];

					for i in range( 0, n ) :
						scalars.SetValue( i, coeffs[ i ] )

				bspline = vtk.vtkBSplineTransform()
				bspline.SetCoefficientData( coefficients )
				bspline.Update()
				transforms.Concatenate( bspline )

			else :

				print( "Error : transform type ", ttype, " not supported" )
				exit( 1 )

	return transforms
