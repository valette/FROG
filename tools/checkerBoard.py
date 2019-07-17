import SimpleITK as sitk
import sys

gridSpacing = 30;

file = sys.argv[1]
image = sitk.ReadImage( file )

size = image.GetSize();
black = sitk.Image( size, sitk.sitkUInt8 )
black.SetOrigin( image.GetOrigin() )
spacing = image.GetSpacing();
black.SetSpacing( spacing )
black.SetDirection( image.GetDirection() )

threshold = sitk.ThresholdImageFilter()
threshold.SetLower( 10 )
threshold.SetOutsideValue ( 100 )
white = threshold.Execute( black )

threshold.SetOutsideValue ( 50 )
grey = threshold.Execute( black )

checker = sitk.CheckerBoardImageFilter();
pattern = [ 0, 0, 0 ]


for i in [ 0, 1, 2 ] :
  pattern[ i ] = int( size[ i ] * spacing[ i ] / gridSpacing )

pattern[ 0 ] = 1

print pattern
checker.SetCheckerPattern( pattern );
board = checker.Execute( grey, white );
sitk.WriteImage( board , "output.nii.gz" )
