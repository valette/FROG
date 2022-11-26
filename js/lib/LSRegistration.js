"use strict";

//const { require } = self;

function getLSRegistration() {

const async = require ('async');
var _ = require ('lodash');

let desk;
try {
    desk = require ('desk-client');
} catch ( e ) {
    desk = window.desk;
}
var getSolver;

if (typeof process !== "undefined") {
	getSolver = () => require(__dirname + '/laplaceSolver.js');
	//desk.include(__dirname +  '/../../frog.json');
} else {
	getSolver = () => laplaceSolver.worker();
}

var LSRegistration = function () {

this.defaultParameters = {
    use3D : true,

    match3d : {
        RansacDist : 40,
        MatchingScale : 1.5,
        MatchingDist : 0.3,
        MatchingDist2 : 0.98,
        ComputeBoundingBoxes : 1,
        RansacMinInliers : 1 //[5, 20, 25, 30, 40]
    },
    surf3d : {
        threshold : 8000,
        spacing: 1.5
    },
    match2d : {
//        RANSAC_iterations : 2000,
        spacing : 2
//        maximum_dimension :150
    },
    edgeRemovalRatio : 0.04,
    finalEdgesRatio : 3
};


var volumes;
var matches;
var result;
var options;

this.register = function (vols, opts, callback) {
	if (typeof opts === 'function') {
		callback = opts;
		opts = {};
	}
    volumes = vols.map(volume => volume instanceof Object ? volume : {file : volume});
    options = _.merge(this.defaultParameters, opts);
    async.waterfall([
        getMatches,
        computeLSRegistration,
    ], callback);
};

var pointsFile = [];
var startTime;

function surf2d (volume, callback) {
    var index = volumes.indexOf(volume);
    async.eachLimit(volumes, 4, function (vol, callback) {
        var index2 = volumes.indexOf(vol);
        if ( index >= index2) {
            callback();
            return;
        }

        var bounds = {};
        var box = vol.boundingBox;
        if (box) {
            bounds.xmin = box[0];
            bounds.xmax = box[1];
            bounds.ymin = box[2];
            bounds.ymax = box[3];
            bounds.zmin = box[4];
            bounds.zmax = box[5];
        }

        var par = Object.assign({ action : "MATCH2D", 
            input_volume1 : vol.file,
            input_volume2 : volume.file
        }, options.match2d, bounds);

        desk.Actions.execute(par, function (err, res) {
                desk.FileSystem.readFile(res.outputDirectory + "/transform.json",
                function (err, res) {
                    if (err) {
                        callback(err);
                        return;
                    }
					res = JSON.parse(res);
                    if (!matches[index]) {
                        matches[index]= [];
                    }
                    matches[index][index2] = res;
                    res.scale = 1 / res.scale;
                    res.translation = res.translation.map(function (c) {
                        return -c;
                    });
                    var temp = res.bboxA;
                    res.bboxA = res.bboxB;
                    res.bboxB = temp;

                    if (res.scale < 0) {
                        res.fail = true;
                    }
                    callback();
                });
            }
        );
    }, callback);
}

function getMatches ( cb ) {
    matches = volumes.map( () => [] );

	if ( options.use3D ) {

		get3DMatches( cb );

	} else {

		async.eachSeries( volumes, surf2d, function ( err ) {
			cb ( err, matches );
		});

	}
}

function get3DMatches( callback ) {

	async.mapLimit( volumes, 3, function ( volume, callback ) {

        desk.Actions.execute( Object.assign( {

                action : "SURF3D",
                input_volume : volume.file,
                writeJSON : 1,
                writeCSVGZ : 0,
                silent : true

            }, options.surf3d ), function ( err, response ) {

                callback( err, response.outputDirectory + "points.json" );

            }
        );

	}, function ( err, result ) {
	    pointsFile = result;
	    var tasks = [];
	    volumes.forEach( function ( volume, index ) {
            for ( var j = 0; j < index; j++ ) {
                tasks.push( { i : index, j : j } );
            }
	    });

        async.each( tasks, function ( task, callback ) {

            desk.Actions.execute(
                Object.assign( { action : "MATCH3D", 
                    input_volume1 : pointsFile[ task.j ], 
                    input_volume2 : pointsFile[ task.i ],
                    stdout : true,
                    silent : true
                },  options.match3d ),

                function ( err, response ) {
					if ( err ) {
						callback( err );
						return;
					}
					matches[ task.j ][ task.i ] = JSON.parse( response.stderr );
					callback();
                }
            );
        }, function ( err ) {
            callback ( err, matches );
        } );

	});
}

function computeLSRegistration (res, callback) {
    matches = res;
    var solver = getSolver();
    var numberOfEdges;
    var ok;

    solver.setInput(matches, options.edgeRemovalRatio, volumes, next);

    function next() {
        async.doWhilst(
            function (callback) {
                solver.iterate(true, function (state) {
                    if (state.error) {
                        ok = false;
                        console.log("error : " + state.error);
                        callback(state.error);
                        return;
                    }
                    result = state;
                    numberOfEdges = state.numberOfEdges;
                    ok = numberOfEdges > matches.length * options.finalEdgesRatio;
                    callback();
                });
            },
            function () {
                return ok; 
        }, function (err) {
			// convert typed arrays to arrays
			if (result) {
				result.valences = Array.prototype.slice.call(result.valences);
				result.positions = Array.prototype.slice.call(result.positions);
			}
            callback (err, result);
        });
    }
}

};

return LSRegistration;
}

const LSRegistration = getLSRegistration();

if ( typeof define === 'function' && define.amd ) {

		define( 'LSRegistration', LSRegistration );

} else if ( 'undefined' !== typeof exports && 'undefined' !== typeof module ) {

		module.exports = LSRegistration;

} else {
    self.LSRegistration = LSRegistration;
}

