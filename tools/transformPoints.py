import argparse
from transformIO import readTransform
import os
from os import listdir
from os.path import isdir, join

parser = argparse.ArgumentParser( description = "Transforms points" , formatter_class=argparse.ArgumentDefaultsHelpFormatter )
parser.add_argument( 'inputPoints', help = 'input points' )
parser.add_argument( 'transform', help = 'transform to apply' )
args = parser.parse_args()

transform = readTransform( args.transform )
file = open( args.inputPoints, mode = 'r')
outputPoints = [];

for line in file.readlines() :
	coords = [ float( n.strip() ) for n in line.split( "," ) ]
	coords2 = transform.TransformPoint( coords )
	output = str( coords2[ 0 ] ) + "," + str( coords2[ 1 ] ) + "," + str( coords2[ 2 ] )
	outputPoints.append( output )

with open( 'output.xyz', 'w') as f:
    f.write( '\n'.join( outputPoints ) )
