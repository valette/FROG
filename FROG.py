import argparse
import nibabel as nib
import os
from os import listdir
from os.path import abspath, dirname, isdir, join, normpath
import time

startTime = time.time()
cwd = os.getcwd()
volumesList = "volumes.txt"
pairsFile = "pairs.bin"
frogPath = normpath( join( dirname( __file__ ), "bin" ) )

parser = argparse.ArgumentParser( description = 'Register a group of volumes', formatter_class=argparse.ArgumentDefaultsHelpFormatter )
parser.add_argument( '-a', dest = 'averageImageSpacing', type = float, help = 'spacing for average image. Not computed if not specified')
parser.add_argument( '-ao', dest = 'averageImageOnly', help = 'only compute average image, skip registration',  action="store_true" )
parser.add_argument( '-i', dest = 'input', help = 'input list of files or directory', required = True )
parser.add_argument( '-dl', dest = 'nDeformableLevels', help = 'number of deformable levels for registration', type = int )
parser.add_argument( '-l', dest = 'limit', type = int, help = 'limit to number of input files' )
parser.add_argument( '-landmarks', dest = 'landmarks', help = 'path to landmarks file' )
parser.add_argument( '-o', dest = 'outputDirectory', help = 'outputDirectory' )
parser.add_argument( '-p', dest = 'numberOfPoints', type = int, help = 'number of keypoints to extract with SURF3D', default = 20000 )
parser.add_argument( '-s', dest = 'spacing', type = float, help = 'spacing for SURF3D', default = 0.75 )
parser.add_argument( '-t', dest = 'threshold', type = float, help = 'detector threshold for SURF3D', default = 0 )
parser.add_argument( '-ras', dest = 'flipToRAS', help = 'ensure RAS orientation for input images',  action="store_true" )

args = parser.parse_args()

def separate():
	print( "******************************************************" )

def execute( cmd ):
	start_time = time.time()
	separate()
	print( "Executing : " + cmd )
	code = os.system( cmd )
	print( "Command : " + cmd )
	print( "Executed in " + str( round( time.time() - start_time ) ) + "s" )
	print( "Exit code : " + str( code ) )
	if code : raise( OSError( code ) )

def flipAndSaveToRAS( filename ):
    
    #Recover the image object
    imageObj = nib.load( filename )
    
    #Get the current orientation
    CurrentOrientation = nib.aff2axcodes(imageObj.affine)
    print("The current orientation is : ", CurrentOrientation)
    
    #Check if the current orientation is already RAS+
    if CurrentOrientation == ('R', 'A', 'S') :
        
        print("Image already recorded into the RAS+ orientation, nothing to do")
        return filename
        
    else :
        #Flip the image to RAS
        flippedImage = nib.as_closest_canonical(imageObj)
                
        ##Check the new orientation
        NewOrientation = nib.aff2axcodes(flippedImage.affine)
        
        #Set Qcode to 1 that the Qform matrix can be used into the further processing
        flippedImage.header['qform_code'] = 1
        
        #Save the flipped image
        nib.save(flippedImage, "RAS.nii.gz")
        
        print("The new orientation is now : ", NewOrientation)
        return join( cwd, "RAS.nii.gz" )

def computeAverageImage( images ) :
	if not args.averageImageSpacing : return
	separate()
	print( "Compute average image" )
	startTime = time.time()
	dummyBin = join( frogPath, "DummyVolumeGenerator" )
	execute( " ".join( [ dummyBin, "bbox.json", str( args.averageImageSpacing ) ] ) )
	transformedImages = []

	for i, image in enumerate( images ) :
		transformBin = join( frogPath, "VolumeTransform" )
		transformedImage = "transformed" + str( i ) + ".nii.gz"
		execute( " ".join( [ transformBin, image, "dummy.mhd", "-t transforms/" + str( i ) + ".json -o " + transformedImage ] ) )
		transformedImages.append( transformedImage )

	averageBin = join( frogPath, "AverageVolumes" )
	execute( " ".join( [ averageBin, " ".join( transformedImages ) ] ) )
	print( "Average image computed in " + str( round( time.time() - startTime ) ) + "s" )


files = []

if isdir( args.input ) :
	for f in sorted( listdir( args.input ) ) :
		for ext in [ ".nii.gz", ".mhd", ".csv.gz" ] :
			if f.endswith( ext ) : files.append( abspath( join( args.input, f ) ) )
else:
	f = open( args.input,mode = 'r' )
	for element in f.read().split( "\n" ) :
		for ext in [ ".nii.gz", ".mhd", ".csv.gz" ] :
			if element.endswith( ext ) : files.append( element )
		
	f.close()

if ( args.limit ) : files = files[ :args.limit ]
print( "There are " + str( len( files ) ) + " files to register : " )
for f in files : print( f )


if args.outputDirectory :
	if not os.path.exists( args.outputDirectory ):
		print( "Output directory does not exist, create it : " + args.outputDirectory )
		os.mkdir( args.outputDirectory )
	os.chdir( args.outputDirectory )

if args.averageImageOnly:
	computeAverageImage( files )
	exit( 0 )

#### compute input volume keypoints if needed
keypointFiles = []

for f in files :
	if f.endswith( ".csv.gz" ) :
		keypointFiles.append( f )
		continue

	separate()
	print ( "Extracting points from " + f )
	if args.flipToRAS: f = flipAndSaveToRAS( f )
	pointsFile = "points" + str( len ( keypointFiles ) )
	surfBin = join( frogPath, "surf3d" )
	execute( surfBin + " " + f + " -s " + str( args.spacing ) + " -t " + str( args.threshold ) + " -n " + str( args.numberOfPoints ) + " -o " + pointsFile )
	keypointFiles.append( join( os.getcwd(),  pointsFile + ".csv.gz") )

separate()

#### compute pairs
volumes = open( volumesList, "w" )
volumes.write( "\n".join( keypointFiles ) )
volumes.close()
matchBin = join( frogPath, "match" )
matchCmd = " ".join( [ matchBin, volumesList, "-o", pairsFile, "-d 1" ] )
execute( matchCmd )

#### register
frogBin = join( frogPath, "frog" )
frogArgs = [ frogBin, "pairs.bin" ]
if args.nDeformableLevels : frogArgs.extend( [ "-dl", str( args.nDeformableLevels ) ] )
if args.landmarks : frogArgs.extend( [ "-l", args.landmarks ] )
execute( " ".join( frogArgs ) )

separate()
print( "Registration done in " + str( round( time.time() - startTime ) ) + "s" )

computeAverageImage( files )
