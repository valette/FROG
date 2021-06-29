import argparse
import json

parser = argparse.ArgumentParser( description = 'Trims a tranfsorm : takes only the n first levels', formatter_class=argparse.ArgumentDefaultsHelpFormatter )
parser.add_argument( 'input', help = 'input transform' )
parser.add_argument( 'n', help = 'number of levels to keep', type = int )
args = parser.parse_args()

with open( args.input ) as f:
	data = json.load( f )

with open('output.json', 'w') as output:
	json.dump( { "transforms" : data[ "transforms" ][ : args.n ] }, output )