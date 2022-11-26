"use strict";

( function () {
    
const { THREE, _, async, desk, EventEmitter } = window;
const FROG = {};

FROG.computeRigidGroupwiseRegistration = async function ( files, options ) {

    const registration = await desk.Actions.executeAsync( {

        action : 'LSRegistration',
        files : await desk.FileSystem.writeCachedFileAsync( 'volumes.json', JSON.stringify( files ) )

    } );

    const result = await desk.FileSystem.readFileAsync( registration.outputDirectory + "registration.json" );
    const positions = JSON.parse( result ).positions;

    return files.map( ( file, index ) => ( {

        volume : file,
        translation : positions.slice( 4 * index, 4 + 4 * index )

    } ) );

};

FROG.DeformableGroupwiseRegistration = function ( volumes, options ) {

    EventEmitter.call( this );
    Object.assign( this, EventEmitter.prototype);

    this.volumes = volumes;
    this.options = options || {};

};

FROG.DeformableGroupwiseRegistration.prototype.execute = async function () {

    let volumes = this.volumes;
    const options = this.options;

    const globalParams = options.globalParams || {};
    const RAWParams = options.RAWParams || {};
    const SURF3DParams = options.SURF3DParams || {};
    const SIFT3DParams = options.SIFT3DParams || {};
    const matchParams = options.matchParams || {};
    const registrationParams = options.registrationParams || {};

    const useSURF = ( options.useSURF !== undefined ) ? options.useSURF : true;
    const useSIFT = ( options.useSIFT !== undefined ) ? options.useSIFT : false;
    const useRAW = ( options.useRAW !== undefined ) ? options.useRAW : false;

    const rootDir = ( await desk.Actions.executeAsync( {

        action : "getRootDir",
        stdout : true

    } ) ).stdout;

    const actions = [];

    if ( useRAW ) {

        actions.push( Object.assign( { action : "SURF3D" },  RAWParams,
            { type : 1 } ) );

    }

    if ( useSURF ) {

        actions.push( Object.assign( { action : "SURF3D" },  SURF3DParams ) );

    }

    if ( useSIFT ) {

        actions.push( Object.assign( { action : "SIFT3D" },  SIFT3DParams ) );

    }

    if ( options.customPoints ) {

        actions.push( { action : "custom", customPoints : options.customPoints } );

    }

    let matchs = [];
    const pointSets = [];

    for ( let action of actions ) {

        let count = 0;

        const points = await async.mapLimit( volumes, 3, async volume => {

            const extractor = action.customPoints ? {

                    outputPoints : action.customPoints[ volumes.indexOf( volume ) ]

                } : await desk.Actions.executeAsync( Object.assign( {
                    inputVolume : volume.volume }, action, globalParams ) );

            let outputPoints = extractor.outputPoints;

            const output = {

                points : extractor.outputPoints,
                extractor : extractor,
                csv : rootDir + outputPoints + ',' + volume.translation.join( ',' )

            };

            if ( action.action === "SURF3D" && options.useTFModel ) {

                const descriptor = await desk.Actions.executeAsync( {
                    action : "applyTFModel", inputCSV : extractor.outputPoints
                } );

                outputPoints = descriptor.outputPoints;
                output.descriptor = descriptor;

            }

            count++;

            this.emit( 'log', action.action + ' keypoints: ' +
                count + '/' + volumes.length + ' done' );

            return output;

        } );

        pointSets.push( points );

        const csvFile = await desk.FileSystem.writeCachedFileAsync(
            "volumes.csv", points.map( obj => obj.csv ).join( "\n" ) );

        this.emit( 'log', 'Computing ' + action.action + ' matches' );

        const match = await desk.Actions.executeAsync( Object.assign( {

            action : "match",
            inputFile: csvFile,
            stdout : true,
            cache : points.map( points => points.points )

        }, matchParams, globalParams ), { listener : mes => this.emit( 'matchLog', mes.data ) }

        );

        this.emit( 'log', 'Total pairing time (ms): ' + match.duration );
        this.emit( 'log' , match.stdout.split( '\n').filter( line => line.includes( 'Nb Match' ) )[ 0 ] );

        matchs.push( match );
    }

    let match = matchs[ 0 ];
    let points = pointSets[ 0 ];

    if ( matchs.length > 1 ) {

        match = await desk.Actions.executeAsync( {

            action : "mergePairs",
            pairs1 : matchs[ 0 ].outputPairs,
            pairs2 : matchs[ 1 ].outputPairs

        } );

    }

    this.emit( 'log', 'Matches computed' );
    this.emit( 'log', 'Computing deformable registration' );

    const registration = await desk.Actions.executeAsync( Object.assign( {

        action : "frog",
        inputPairs : match.outputPairs,
        nImages : volumes.length,

    }, registrationParams, globalParams ), { listener : mes => this.emit( 'registrationLog', mes.data ) } );

    this.emit( 'log', 'Registration done' );
    this.emit( 'log', 'Total registration time (ms): ' + registration.duration );
    const output = _.cloneDeep( volumes );
    output.match = match;

    output.forEach( ( file, index ) => {

        file.points = points[ index ].points;
        file.extractor = points[ index ].extractor;
        file.transform = registration.outputDirectory + "transforms/" + index + '.json';

    } );

    await desk.FileSystem.writeJSONAsync( registration.outputDirectory +
        "registration.json", output );

    return { match, registration, volumes : output };
};


FROG.CommonSpaceMeanImage = function ( opts ) {

    EventEmitter.call( this );
    Object.assign( this, EventEmitter.prototype);
    this.opts = opts || {};

};

FROG.CommonSpaceMeanImage.prototype.execute = async function () {

    const opts = this.opts;
    const registration = opts.registration;

    const bbox = JSON.parse( await desk.FileSystem.readFileAsync(
        registration.outputDirectory + "bbox.json" ) ).bbox;

    let min = new THREE.Vector3().fromArray( bbox[ 0 ] );
    let max = new THREE.Vector3().fromArray( bbox[ 1 ] );

    let target;

    if ( opts.targetAverageVolume ) {

        target =  opts.targetAverageVolume;
        let samplingRatio = opts.samplingRatio || 1;

        if  ( samplingRatio < 1 ) {

            target = ( await desk.Actions.executeAsync( {

                action : "volume_subsample",
                input_volume : target,
                ratio : samplingRatio

            } ) ).outputDirectory + "out.mhd";

        }

    } else {

        const spacing = opts.spacing || 2;
        max.sub( min ).divideScalar( spacing ).ceil();

        const dummyVolumeGenerator = await desk.Actions.executeAsync( {

            action : "DummyVolumeGenerator",
            ox : min.x,
            oy : min.y,
            oz : min.z,
            sx : spacing,
            sy : spacing,
            sz : spacing,
            dx : max.x,
            dy : max.y,
            dz : max.z

        } );

        target = dummyVolumeGenerator.outputDirectory + "volume.mhd";

    }

    let count = 0;
    this.emit( 'log',  'Transforming images. ' );

    const files = opts.volumes.map( volume => volume.volume );
    const transformedVolumes = await async.mapLimit( files, 4, async file => {

        const index = files.indexOf( file );

        const transformVolume = await desk.Actions.executeAsync( {

            action : 'VolumeTransform',
            source : file,
            transform : registration.outputDirectory + "transforms/" + index + '.json',
            reference : target
//            outputFileName : "output.mhd"

        } );

        count++;
        this.emit( 'log', 'Transforming images. ' + count + '/' + files.length + ' done' );
        return transformVolume.outputDirectory + "output.mhd";

    } );


    const averageVolumes = await desk.Actions.executeAsync( {

        action : 'averageVolumes',
        inputVolumes : transformedVolumes,
        normalize : 0

    } );

    return { transformedVolumes,
        averageVolume : averageVolumes.outputDirectory + 'average.nii.gz' };

};

window.FROG = FROG;
} ) ();