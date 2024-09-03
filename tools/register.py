import argparse
from flip_tool import flip_image
import nibabel as nib
import os
from os import listdir
from os.path import abspath, dirname, join, normpath
import time

startTime = time.time()
cwd = os.getcwd()
volumesList = "volumes.txt"
frogPath = normpath( join( dirname( __file__ ), "../bin" ) )

parser = argparse.ArgumentParser( description = 'Register a volume against an already registered group', formatter_class=argparse.ArgumentDefaultsHelpFormatter )
parser.add_argument( '-a', dest = 'all', help = 'register all', action="store_true" )
parser.add_argument( '-d', dest = 'inputDir', help = 'input registered group directory', required = True )
parser.add_argument( '-dl', dest = 'nDeformableLevels', help = 'number of deformable levels for registration', type = int )
parser.add_argument( '-l', dest = 'landmarks', help = 'path to landmarks file' )
parser.add_argument( '-li', dest = 'nLinearIterations', help = 'number of linear iterations for registration', type = int )
parser.add_argument( '-i', dest = 'inputVolume', help = 'input volume', required = True )
parser.add_argument( '--orientation', help = 'Force image orientation' )
parser.add_argument('--fake_orientation', help = "Fake image orientation" )
parser.add_argument( '-p', dest = 'numberOfPoints', type = int, help = 'number of keypoints to extract with SURF3D', default = 20000 )
parser.add_argument( '-n', dest = 'numberOfReferences', type = int, help = 'number of reference volumes to register against' )
parser.add_argument( '-s', dest = 'spacing', type = float, help = 'spacing for SURF3D', default = 0.75 )
parser.add_argument( '-t', dest = 'threshold', type = float, help = 'detector threshold for SURF3D', default = 0 )
parser.add_argument( '-vt', dest = 'volumeThreshold', type = float, help = 'volume Threshold' )
parser.add_argument( '-cmin', type = float, help = 'min clamp image values' )
parser.add_argument( '-cmax', type = float, help = 'max clamp image values' )
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

points = []

def getId( elem ):
	return int( elem[ 6: -7] )

### get reference points
for f in os.listdir( args.inputDir ):
	if not f.endswith( ".csv.gz" ) : continue
	if not f.startswith( "points" ) : continue
	points.append( f )

points.sort( key = getId )
if args.numberOfReferences : points = points[ 0 : args.numberOfReferences ]
print( str( len( points ) ) + " references" )

#### compute input volume keypoints if needed
inputPoints = "";

if not args.inputVolume.endswith( '.csv.gz' ) :
	if args.orientation or args.fake_orientation:
		inputVolume = flip_image( args.inputVolume, threshold = args.volumeThreshold, orientation = args.orientation, fake_orientation = args.fake_orientation )
	else : inputVolume = args.inputVolume
	surfBin = join( frogPath, "surf3d" )
	surfArgs = [ surfBin, inputVolume, "-s", str( args.spacing ), "-t", str( args.threshold ), "-n", str( args.numberOfPoints )]
	if args.cmin != None : surfArgs.extend( [ "-cmin", str( args.cmin ) ] )
	if args.cmax != None : surfArgs.extend( [ "-cmax", str( args.cmax ) ] )
	execute( " ".join( surfArgs ) )
	inputPoints = join( cwd, "points.csv.gz" )
else:
	print ( "Input is a .csv.gz file, no need for extraction" )
	inputPoints = abspath( args.inputVolume )

#### compute pairs
volumes = open( volumesList, "w" )
for pt in points : volumes.write( join( args.inputDir, pt ) + "\n" )
volumes.write( inputPoints )
volumes.close()
matchBin = join( frogPath, "match" )
matchCmd = matchBin + " " + volumesList + " -o pairs.bin -d 1"
if not args.all : matchCmd += " -targ " + str( len( points ) )
execute( matchCmd )

#### register
frogBin = join( frogPath, "frog" )
frogArgs = [ frogBin, "pairs.bin", "-j -fd", args.inputDir + "/transforms" ]
if args.landmarks : frogArgs.extend( [ "-l", args.landmarks ] )
if args.nDeformableLevels != None : frogArgs.extend( [ "-dl", str( args.nDeformableLevels ) ] )
if args.nLinearIterations : frogArgs.extend( [ "-li", str( args.nLinearIterations ) ] )
if not args.all : frogArgs.extend( [ "-fi", str( len( points ) ) ] )
execute( " ".join( frogArgs ) )

separate()
print( "Registration done in " + str( round( time.time() - startTime ) ) + "s" )
separate()
