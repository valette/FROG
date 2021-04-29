import argparse
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
parser.add_argument( '-dl', dest = 'nDeformableLevels', help = 'number of deformable levels for registration', type = int, default = 3 )
parser.add_argument( '-i', dest = 'inputVolume', help = 'input volume', required = True )
parser.add_argument( '-p', dest = 'numberOfPoints', type = int, help = 'number of keypoints to extract with SURF3D', default = 20000 )
parser.add_argument( '-n', dest = 'numberOfReferences', type = int, help = 'number of reference volumes to register against' )
parser.add_argument( '-s', dest = 'spacing', type = float, help = 'spacing for SURF3D', default = 0.75 )
parser.add_argument( '-t', dest = 'threshold', type = float, help = 'detector threshold for SURF3D', default = 0 )
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
	inputVolume = flipAndSaveToRAS( args.inputVolume )
	surfBin = join( frogPath, "surf3d" )
	execute( surfBin + " " + inputVolume + " -s " + str( args.spacing ) + " -t " + str( args.threshold ) + " -n " + str( args.numberOfPoints ) )
	inputPoints = join( cwd, "points.csv.gz" )
else:
	print ( "Input is a .csv.gz file, no need for extraction" )
	inputPoints = abspath( args.inputVolume )

#### compute pairs
volumes = open( volumesList, "w" )
for pt in points:
	volumes.write( join( args.inputDir, pt ) + "\n" )

volumes.write( inputPoints )
volumes.close()

matchBin = join( frogPath, "match" )
matchCmd = matchBin + " " + volumesList + " -o pairs.bin -d 1"
if not args.all : matchCmd += " -targ " + str( len( points ) )
execute( matchCmd )

#### register
frogBin = join( frogPath, "frog" )
frogCmd = frogBin + " pairs.bin -fd " + args.inputDir + "/transforms" + " -dl " + str( args.nDeformableLevels )
if not args.all : frogCmd += " -fi " + str( len( points ) )
execute( frogCmd )

separate()
print( "Registration done in " + str( round( time.time() - startTime ) ) + "s" )
separate()