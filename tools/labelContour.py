import SimpleITK as sitk

import sys

file = sys.argv[1]
image = sitk.ReadImage( file )

contour = sitk.LabelContourImageFilter()
output = contour.Execute( image )

sitk.WriteImage( output , "output.nii.gz" )
sitk.WriteImage( output , "output.mhd" )
