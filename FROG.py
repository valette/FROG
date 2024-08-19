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
regParse.add_argument( '-se', '--skip-existing', dest="skipExisting", help = 'Do not recompute files if they already exist', action = "store_true" )
regParse.add_argument( '-limit', dest = 'limit', type = int, help = 'limit number of input files' )
regParse.add_argument( '-o', dest = 'outputDirectory', help = 'outputDirectory' )
regParse.add_argument( '-j', dest = 'useSingleJSONTransformFile', help = 'use a single JSON file to store each transform', action="store_true"  )
frogParse = parser.add_argument_group('Registration options')
frogParse.add_argument( '-dl', dest = 'deformableLevels', help = 'number of deformable levels', type = int )
frogParse.add_argument( '-di', dest = 'deformableIterations', help = 'number of deformable iterations per level', type = int )
frogParse.add_argument( '-g', dest = 'gridSpacing', type = float, help = 'initial grid spacing' )
frogParse.add_argument( '-lanchor', dest = 'linearInitializationAnchor', help = 'Linear initialization Anchor', type = float, nargs = 3, default = [0.5,0.5,0.5] )
frogParse.add_argument( '-li', dest = 'linearIterations', help = 'number of linear iterations', type = int )
frogParse.add_argument( '-ri', dest = 'RANSACIterations', help = 'number of RANSAC iterations', type = int )
frogParse.add_argument( '-l', dest = 'landmarks', help = 'path to landmarks file' )
frogParse.add_argument( '-il', dest = 'invertLandmarks', help = 'revert landmarks coordinates', type = int, default = 1 )
frogParse.add_argument( '-wp', dest = 'writePairs', help = 'write list of pairs to file', action="store_true" )
matchParser = parser.add_argument_group('Match options')
matchParser.add_argument( '-md', dest = 'matchDistance', type = float, default = 10000000000, help = 'maximum descriptor distance' )
matchParser.add_argument( '-d2', dest = 'ratio', type = float, default = 1, help = 'maximum second match distance ratio' )

SURFParser = parser.add_argument_group('SURF3D options')
SURFParser.add_argument( '-m', dest = 'masks', help = 'path to masks' )
SURFParser.add_argument( '-cmin', type = float, help = 'min clamp image values' )
SURFParser.add_argument( '-cmax', type = float, help = 'max clamp image values' )
SURFParser.add_argument( '-p', dest = 'numberOfPoints', type = int, help = 'number of keypoints to extract', default = 20000 )
SURFParser.add_argument( '-pad', dest = 'padding', type = float, help = 'mirror padding length' )
SURFParser.add_argument( '-s', dest = 'spacing', type = float, help = 'spacing', default = 0.75 )
SURFParser.add_argument( '-t', dest = 'threshold', type = float, help = 'detector threshold', default = 0 )
SURFParser.add_argument( '-ras', dest = 'flipToRAS', help = 'ensure RAS orientation for input images',  action="store_true" )
SURFParser.add_argument( '-rast', dest = 'useTempd', help = 'use in-memory temporary directory when converting to RAS',  action="store_true" )
averageParse = parser.add_argument_group('Average image computing')
averageParse.add_argument( '-a', dest = 'imageSpacing', type = float, help = 'spacing for average image. Not computed if not specified')
averageParse.add_argument( '-ao', dest = 'averageImageOnly', help = 'only compute average image, skip registration',  action="store_true" )
averageParse.add_argument( '-ad', dest = 'averageImageDirectory', help = 'Average image subdirectory' )
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
        img_data = flippedImage.get_fdata()
        img_conv = nib.Nifti1Image(img_data.astype(flippedImage.header.get_data_dtype()), flippedImage.affine, flippedImage.header)

        #Set Qcode to 1 that the Qform matrix can be used into the further processing
        img_conv.header['qform_code'] = 1
        
        #Save the flipped image
        nib.save(img_conv, RASFile )

        print("The new orientation is now : ", NewOrientation)
        return RASFile

def computeAverageImage( images ) :
	if not args.imageSpacing : return
	separate()
	print( "Compute average image" )
	startTime = time.time()
	dummyBin = join( frogPath, "DummyVolumeGenerator" )
	execute( " ".join( [ dummyBin, "bbox.json", str( args.imageSpacing ) ] ) )
	dummyFile = "dummy.mhd";
	if args.averageImageDirectory:
		if not os.path.exists( args.averageImageDirectory ):
			os.mkdir( args.averageImageDirectory )
		for file in [ "dummy.mhd", "dummy.zraw" ]:
			os.rename( file, join( args.averageImageDirectory, file ) )
		dummyFile = join( args.averageImageDirectory, dummyFile )
	transformedImages = []

	for i, image in enumerate( images ) :
		separate()
		transformBin = join( frogPath, "VolumeTransform" )
		transformedImage = "transformed" + str( i ) + ".nii.gz"
		if args.averageImageDirectory:
			transformedImage = join( args.averageImageDirectory, transformedImage )
		if args.flipToRAS: image = flipAndSaveToRAS( image )
		execute( " ".join( [ transformBin, image, dummyFile, "-t transforms/" + str( i ) + ".json -o " + transformedImage ] ) )
		transformedImages.append( transformedImage )

	averageBin = join( frogPath, "AverageVolumes" )
	execute( " ".join( [ averageBin, " ".join( transformedImages ) ] ) )
	print( "Average image computed in " + str( round( time.time() - startTime ) ) + "s" )
	if args.averageImageDirectory:
		for file in [ "average.nii.gz", "stdev.nii.gz" ]:
			os.rename( file, join( args.averageImageDirectory, file ) )

def getFileList( inputPath ) :

	files = []

	if isdir( inputPath ) :
		for f in sorted( listdir( inputPath ) ) :
			for ext in [ ".nii.gz", ".mhd", ".csv.gz" ] :
				if f.endswith( ext ) : files.append( abspath( join( inputPath, f ) ) )
	else:
		f = open( inputPath, mode = 'r' )
		for element in f.read().split( "\n" ) :
			if element.startswith( "#" ) : continue
			for ext in [ ".nii.gz", ".mhd", ".csv.gz" ] :
				if element.endswith( ext ) : files.append( join( dirname( inputPath) , element ) )

		f.close()

	return files


files = getFileList( args.input )
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
if args.masks: maskFiles = getFileList( args.masks )


#### compute input volume keypoints if needed
keypointFiles = []

for index, f in enumerate( files ):
	if f.endswith( ".csv.gz" ) :
		keypointFiles.append( f )
		continue

	separate()
	pointsFile = "points" + str( len ( keypointFiles ) )
	fullPointsFile = join( os.getcwd(),  pointsFile + ".csv.gz")
	keypointFiles.append( fullPointsFile )
	if ( args.skipExisting and os.path.exists( fullPointsFile ) ):
		print( "Points file ", fullPointsFile, "already exists, skipping" )
		continue
	print ( "Extracting points from " + f )
	if args.flipToRAS: f = flipAndSaveToRAS( f )
	surfBin = join( frogPath, "surf3d" )
	surfArgs = [ surfBin, f, "-s", str( args.spacing ), "-t", str( args.threshold ), "-n", str( args.numberOfPoints ), "-o", pointsFile ]
	if len( maskFiles ) : surfArgs.extend( [ "-m", maskFiles[ index ] ] )
	if args.cmin != None : surfArgs.extend( [ "-cmin", str( args.cmin ) ] )
	if args.cmax != None : surfArgs.extend( [ "-cmax", str( args.cmax ) ] )
	if args.padding != None : surfArgs.extend( [ "-pad", str( args.padding ) ] )
	execute( " ".join( surfArgs ) )

separate()

#### compute pairs
if args.skipExisting and os.path.exists( pairsFile ):
	print( "Pairs file pairs.bin already exists, skipping computation" )
else:
	volumes = open( volumesList, "w" )
	volumes.write( "\n".join( keypointFiles ) )
	volumes.close()
	matchBin = join( frogPath, "match" )
	matchCmd = " ".join( [ matchBin, volumesList, "-o", pairsFile, "-d", str( args.matchDistance ), "-np", str( args.numberOfPoints ), "-d2", str( args.ratio) ] )
	execute( matchCmd )

#### register
frogBin = join( frogPath, "frog" )
frogArgs = [ frogBin, "pairs.bin" ]
if args.deformableLevels != None : frogArgs.extend( [ "-dl", str( args.deformableLevels ) ] )
if args.RANSACIterations : frogArgs.extend( [ "-ri", str( args.RANSACIterations ) ] )
if args.linearIterations : frogArgs.extend( [ "-li", str( args.linearIterations ) ] )
if args.deformableIterations : frogArgs.extend( [ "-di", str( args.deformableIterations ) ] )
if args.writePairs : frogArgs.append( "-wp 1" )
if args.landmarks : frogArgs.extend( [ "-l", args.landmarks ] )
if args.gridSpacing : frogArgs.extend( [ "-g", str( args.gridSpacing ) ] )
if args.invertLandmarks : frogArgs.extend( [ "-il", str( args.invertLandmarks ) ] )
if args.useSingleJSONTransformFile : frogArgs.extend( [ "-j" ] )
if args.linearInitializationAnchor : frogArgs.extend( [ "-lanchor", *map( lambda x : str( x ), args.linearInitializationAnchor ) ] )
execute( " ".join( frogArgs ) )

separate()
print( "Registration done in " + str( round( time.time() - startTime ) ) + "s" )

computeAverageImage( files )
