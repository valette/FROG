import SimpleITK as sitk

import sys
finalSpacing = 0.5;
gridSpacing = 30 / finalSpacing;
sigma = 0.5

file = sys.argv[1]
image = sitk.ReadImage( file )
gridSource = sitk.GridImageSource()

gridSource.SetGridSpacing( [ gridSpacing, gridSpacing, gridSpacing ] )
gridSource.SetOrigin( image.GetOrigin() );
spacing = image.GetSpacing();
gridSource.SetSigma( [ sigma, sigma, sigma] );
gridSource.SetDirection( image.GetDirection() );
size = image.GetSize();

s =[ 0,0,0];
sp =[ 0,0,0];

for i in [ 0, 1, 2 ]:
    s[ i ] = int ( round( size[ i ] * spacing[ i ] / finalSpacing ) + 1 );
    sp[ i ] = finalSpacing;

gridSource.SetSize( s );
gridSource.SetSpacing( sp );
output = gridSource.Execute()
sitk.WriteImage( output , "output.nii.gz" )
