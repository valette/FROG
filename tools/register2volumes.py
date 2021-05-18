import argparse
import os
from os import listdir
from os.path import abspath, dirname, isdir, join, normpath
import tempfile
import time

startTime = time.time()
cwd = os.getcwd()
volumesList = "toRegister.txt"
masksList = "masks.txt"
rootPath = normpath( join( dirname( __file__ ), ".." ) )

parser = argparse.ArgumentParser( description = 'Register 2 images', formatter_class=argparse.ArgumentDefaultsHelpFormatter )
regParse = parser.add_argument_group('General options')
regParse.add_argument( 'source', help = 'source image' )
regParse.add_argument( 'target', help = 'target image' )
regParse.add_argument( '-m1', dest = 'mask1', help = 'source mask' )
regParse.add_argument( '-m2', dest = 'mask2', help = 'target mask' )
regParse.add_argument( '-o', dest = 'output', help = 'output image' )
frogParse = parser.add_argument_group('Registration options')
frogParse.add_argument( '-dl', dest = 'deformableLevels', help = 'number of deformable levels', type = int )
frogParse.add_argument( '-di', dest = 'deformableIterations', help = 'number of deformable iterations per level', type = int )
frogParse.add_argument( '-g', dest = 'gridSpacing', type = float, help = 'initial grid spacing' )
frogParse.add_argument( '-li', dest = 'linearIterations', help = 'number of linear iterations', type = int )
frogParse.add_argument( '-l', dest = 'landmarks', help = 'path to landmarks file' )
frogParse.add_argument( '-wp', dest = 'writePairs', help = 'write list of pairs to file', action="store_true" )
SURFParser = parser.add_argument_group('SURF3D options')
SURFParser.add_argument( '-m', dest = 'masks', help = 'path to masks' )
SURFParser.add_argument( '-p', dest = 'numberOfPoints', type = int, help = 'number of keypoints to extract', default = 20000 )
SURFParser.add_argument( '-s', dest = 'spacing', type = float, help = 'spacing', default = 0.75 )
SURFParser.add_argument( '-t', dest = 'threshold', type = float, help = 'detector threshold', default = 0 )
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


def computeAverageImage( images ) :
	if not args.imageSpacing : return
	separate()
	print( "Compute average image" )
	startTime = time.time()
	dummyBin = join( frogPath, "DummyVolumeGenerator" )
	execute( " ".join( [ dummyBin, "bbox.json", str( args.imageSpacing ) ] ) )
	transformedImages = []

	for i, image in enumerate( images ) :
		separate()
		transformBin = join( frogPath, "VolumeTransform" )
		transformedImage = "transformed" + str( i ) + ".nii.gz"
		if args.flipToRAS: image = flipAndSaveToRAS( image )
		execute( " ".join( [ transformBin, image, "dummy.mhd", "-t transforms/" + str( i ) + ".json -o " + transformedImage ] ) )
		transformedImages.append( transformedImage )

	averageBin = join( frogPath, "AverageVolumes" )
	execute( " ".join( [ averageBin, " ".join( transformedImages ) ] ) )
	print( "Average image computed in " + str( round( time.time() - startTime ) ) + "s" )

files = [ abspath( args.source ), abspath( args.target ) ]
masks = []
if args.mask1 and args.mask2 :
	masks = [ abspath( args.mask1 ), abspath( args.mask2 ) ]

print( "There are " + str( len( files ) ) + " files to register : " )
for f in files : print( f )
separate()

frogBin = join( rootPath, "FROG.py" )
volumes = open( volumesList, "w" )
volumes.write( "\n".join( files ) )
volumes.close()

frogArgs = [ "python", frogBin, volumesList, "-s", str( args.spacing ), "-t", str( args.threshold ), "-p", str( args.numberOfPoints ) ]

if len( masks ) :
	print( "Using masks : " )
	for f in masks : print( f )
	mFile = open( masksList, "w" )
	mFile.write( "\n".join( masks ) )
	mFile.close()
	frogArgs.extend( [ "-m", abspath( masksList ) ] )

if args.deformableLevels : frogArgs.extend( [ "-dl", str( args.deformableLevels ) ] )
if args.linearIterations : frogArgs.extend( [ "-li", str( args.linearIterations ) ] )
if args.deformableIterations : frogArgs.extend( [ "-di", str( args.deformableIterations ) ] )
if args.writePairs : frogArgs.append( "-wp" )
if args.landmarks : frogArgs.extend( [ "-l", args.landmarks ] )
if args.gridSpacing : frogArgs.extend( [ "-g", str( args.gridSpacing ) ] )
if args.spacing : frogArgs.extend( [ "-s", str( args.spacing ) ] )
if args.numberOfPoints : frogArgs.extend( [ "-p", str( args.numberOfPoints ) ] )
execute( " ".join( frogArgs ) )


separate()
transformBin = join( rootPath, "bin/VolumeTransform" )
transformArgs = [ transformBin, files[ 0 ], files[ 1 ] ]
transformArgs.append( "-t transforms/0.json" )
transformArgs.append( "-ti transforms/1.json" )
if ( args.output ) : transformArgs.extend( [ "-o", args.output ] )
execute( " ".join( transformArgs ) )

separate()
print( "All done in " + str( round( time.time() - startTime ) ) + "s" )
