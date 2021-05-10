import argparse
import nibabel as nib
import os
from os import listdir
from os.path import abspath, dirname, isdir, join, normpath
import tempfile
import time

startTime = time.time()
cwd = os.getcwd()
volumesList = "volumes.txt"
pairsFile = "pairs.bin"
frogPath = normpath( join( dirname( __file__ ), "bin" ) )

parser = argparse.ArgumentParser( description = 'Register a group of volumes', formatter_class=argparse.ArgumentDefaultsHelpFormatter )
regParse = parser.add_argument_group('General options')
regParse.add_argument( 'input', help = 'input list of files or directory' )
regParse.add_argument( '-limit', dest = 'limit', type = int, help = 'limit number of input files' )
regParse.add_argument( '-o', dest = 'outputDirectory', help = 'outputDirectory' )
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
SURFParser.add_argument( '-ras', dest = 'flipToRAS', help = 'ensure RAS orientation for input images',  action="store_true" )
SURFParser.add_argument( '-rast', dest = 'useTempd', help = 'use in-memory temporary directory when converting to RAS',  action="store_true" )
averageParse = parser.add_argument_group('Average image computing')
averageParse.add_argument( '-a', dest = 'imageSpacing', type = float, help = 'spacing for average image. Not computed if not specified')
averageParse.add_argument( '-ao', dest = 'averageImageOnly', help = 'only compute average image, skip registration',  action="store_true" )
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
        nib.save(flippedImage, RASFile )
        
        print("The new orientation is now : ", NewOrientation)
        return RASFile

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

tempD = 0
RASFile = 0
if args.useTempd :
	tempD = tempfile.TemporaryDirectory();
	RASFile = join( tempD.name, "RAS.nii" )
else:
	RASFile = join( os.getcwd(), "RAS.nii.gz" )

if args.averageImageOnly:
	computeAverageImage( files )
	exit( 0 )

maskFiles = []

if args.masks:
	for f in sorted( listdir( args.masks ) ) :
		for ext in [ ".nii.gz", ".mhd", ".csv.gz" ] :
			if f.endswith( ext ) :
				maskFiles.append( abspath( join( args.masks, f ) ) )
				print (f )


#### compute input volume keypoints if needed
keypointFiles = []

for index, f in enumerate( files ):
	if f.endswith( ".csv.gz" ) :
		keypointFiles.append( f )
		continue

	separate()
	print ( "Extracting points from " + f )
	if args.flipToRAS: f = flipAndSaveToRAS( f )
	pointsFile = "points" + str( len ( keypointFiles ) )
	surfBin = join( frogPath, "surf3d" )
	surfArgs = [ surfBin, f, "-s", str( args.spacing ), "-t", str( args.threshold ), "-n", str( args.numberOfPoints ), "-o", pointsFile ]
	if len( maskFiles ) : surfArgs.extend( [ "-m", maskFiles[ index ] ] )
	execute( " ".join( surfArgs ) )
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
if args.deformableLevels : frogArgs.extend( [ "-dl", str( args.deformableLevels ) ] )
if args.linearIterations : frogArgs.extend( [ "-li", str( args.linearIterations ) ] )
if args.deformableIterations : frogArgs.extend( [ "-di", str( args.deformableIterations ) ] )
if args.writePairs : frogArgs.append( "-wp 1" )
if args.landmarks : frogArgs.extend( [ "-l", args.landmarks ] )
if args.gridSpacing : frogArgs.extend( [ "-g", str( args.gridSpacing ) ] )
execute( " ".join( frogArgs ) )

separate()
print( "Registration done in " + str( round( time.time() - startTime ) ) + "s" )

computeAverageImage( files )
