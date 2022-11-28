'use strict';

const { async, THREE, desk, FROG, qx, visceral, _, ACVD } = window;

const filter1 = "all"; // registration batch. Can be "all", "women" or "men"
const filter2 = "all"; // average image computation. Can be "all", "women" or "men" or a volume id

let spacing = 5; // spacing of the output mean image

let placesAll = [ // input directories

    'data/heads/heads.json',
    'data/heads/heads2.json',
    "data/volumes-01_MHD-250.json",
    'big/01_Registrations/registration-001-931.json',
    'big/visceral/volumes',
    'big/visceral/silver/volumes/CTce_ThAb',
    'big/visceral/volumes_CT_wb',
    'big/LPBA40/LPBA40_mri.json',
    'data/loiseau/epauleDroiteGauche_Frog'

];

let displayLandmarks = false;

let useRAW = false;
let useSIFT = false;
let useSURF = true;

let defaultNumberOfKeypoints = 20000;

const RAWParams = {

    spacing : 0.75,
    threshold : 0,

};

const SURF3DParams = {

//forceUpdate : true,
//    spacing : 0.75,
    spacing : 0.75,
    threshold : 0,
//    radius : 5

};


const SIFT3DParams = {

    spacing : 0.75,
    threshold : 0

};

const matchParams = {

    distance : 0.5,
//    all : 1,
    distanceToSecond : 0.998

};

const registrationParams = {
//forceUpdate : true,
displayLinear : 1,
//statsInterval : 1,
//statsInterval : 1,
//    linearIterations : 200,
    deformableLevels : 3,
    deformableIterations : 200,
//    deformableAlpha : 0.3,
//    initialGridSize : 100,
//    inlierThreshold : 0.01

// old matlab stuff
//    nIterations : 799,
//    iterationsPerLevel : 200

};

// meshes to extract from average image
let organs = {

    bones: {

        threshold : 200,

    },

    skin : {

        threshold : -200,
        transparent : true,
        renderOrder : 1000,

        material : {

            side : THREE.DoubleSide,
            opacity : 0.7,
            color : 'rgb(255,192,203)',
            transparent : true

        }

    }

};

////////////////////////////////////////////////////////////////////////
let files = []; // this array will be filled with volume files to register

function isVolume( file ) {

    file = file.toLowerCase();
    return ( file.endsWith( ".mhd" ) || file.endsWith( ".nii.gz" ) );

}

async function run() {

let targetAverageVolume; // leave blank for automatic bounds computation
//targetAverageVolume = 'big/visceral/registrations/Ours20/output.mhd';
//targetAverageVolume = "data/heads/head.mhd";
const samplingRatio = 0.5;

var size = 10000; // maximal global size
//size =3
const spacingX = 500; // grid x spacing for display
const spacingZ = 1500; // grid z spacing for display

const loadVolumeConcurrency = 5;

let width;

// setup UI
const viewer = new desk.MeshViewer( {

    orthographic : true,
    cameraFront : [ 0, 1, 0 ],
    cameraUp : [ 0, 0, 1 ]

} );

let vv;
viewer.getControls().noRotate = true;

//if ( desk.auto ) viewer.fillScreen();

const list = new qx.ui.form.List();
viewer.add( list, { top : 0, left : 0 } );
list.set( { width : 300, selectionMode : 'multi' } );
for ( let place of placesAll ) list.add( new qx.ui.form.ListItem( place ) );
list.setSelection( [ list.getChildren()[ 4 ] ] );
list.setDroppable(true);
list.addListener( "drop", e => {

	if ( e.supportsType( "file" ) ) list.add(
	    new qx.ui.form.ListItem( e.getData( "file" ) ) );

} );


const startButton = new qx.ui.form.Button( 'Start' );
viewer.add( startButton, { top : 0, left : 300 } );
if ( !desk.auto ) await new Promise( resolve => startButton.addListenerOnce( 'execute', resolve ) );

const places = list.getSelection().map( item => item.getLabel() );
startButton.destroy();
list.destroy();

viewer.getRenderer().sortObjects = false;
const statusLabel = new qx.ui.basic.Label( '' );
statusLabel.set( { backgroundColor : "white" } );
statusLabel.setValue = _.throttle( statusLabel.setValue.bind( statusLabel ), 500 );
viewer.add( statusLabel, { left : "15%", bottom : 50 } );

let presets;

// build file list
for ( let place of places ) {

    if ( place.endsWith( ".json" ) ) {

        presets = JSON.parse( await desk.FileSystem.readFileAsync( place ) );
        const list = Array.isArray( presets ) ? presets : presets.volumes;
        files = files.concat( list );

    } else if ( isVolume( place ) ) {

        files.push( place );

    } else {

        await desk.FileSystem.traverseAsync( place, file => {
            if ( isVolume( file ) ) files.push( file ); } );

    }
}

files = _.uniq( files )
    .filter( getVolumeFilter( filter1 ) );

files.length = Math.min( files.length, size );
size = files.length;
width =  files.length;

// detect if we have landmarks for the files

if ( window.visceral ) {

    const file = files[ 0 ];
    const id = visceral.getFileId( file );
    if ( file.includes( "visceral" ) && !isNaN( id ) ) {

        registrationParams.landmarks = desk.FileSystem.getFileDirectory( visceral.getLandmarkFile( id ) );

    }

}

// Load all volumes

let numberOfLoadedVolumes;

const meshes = await async.mapLimit( files, loadVolumeConcurrency,

    async function( file ) {

        const id = files.indexOf( file );
        const anchor = new THREE.Group();

        viewer.addMesh( anchor, { label : '' + id } );

        const mesh = await viewer.addVolumeAsync( file , {

            orientations : [ 2 ],
            sliceWith : { sampling_ratio : 0.1 },
            parent : anchor

        } );

        mesh.userData.id = id;
        mesh.userData.file = file;

        const i = id % width;
        const j = Math.floor( id / width );

        anchor.position.x = i * spacingX;
        anchor.position.z = - j * spacingZ;
        anchor.userData.initial = anchor.position.clone();
        anchor.userData.final = new THREE.Vector3();

        viewer.resetView();
        if ( !numberOfLoadedVolumes ) numberOfLoadedVolumes = 1;
        statusLabel.setValue( numberOfLoadedVolumes++ + " volumes loaded" );

        return mesh;

    }
);

statusLabel.setValue( "All " + meshes.length + " volumes loaded." );

// compute rigid registration ( or read it from file )

let positions;
if ( presets && presets.positions ) {

    positions = presets.positions;

    for ( let mesh of meshes ) {

        mesh.position.fromArray( positions, 4 * presets.volumes.indexOf( mesh.userData.file ) );

    }

} 

viewer.render();
statusLabel.setValue( "Click on 'register' to start registration " );

for ( let mesh of meshes ) {

    const bounds = mesh.children[ 0 ].userData.viewerProperties.volumeSlice.getBounds();
    mesh.userData.zMin = bounds [ 4 ];
    mesh.userData.zMax = bounds [ 5 ];

}

const container = new qx.ui.container.Composite( new qx.ui.layout.VBox() );
container.setDecorator('main');
const container2 = new qx.ui.container.Composite( new qx.ui.layout.VBox() );
container2.setVisibility( 'excluded' );

container.set( {

    width : 200,
    backgroundColor : 'white'

} );

viewer.add( container, { top : 0, right : 10 } );

const label = new qx.ui.basic.Label( '' );
container.add( label );

container.add( container2 );

const useRAWBox = new qx.ui.form.CheckBox( 'use RAW' );
useRAWBox.setValue( useRAW );
useRAWBox.addListener( "changeValue", e => useRAW = e.getData() );
container.add( useRAWBox );

const useSURFBox = new qx.ui.form.CheckBox( 'use SURF' );
useSURFBox.setValue( useSURF );
useSURFBox.addListener( "changeValue", e => useSURF = e.getData() );
container.add( useSURFBox );

const useSIFTBox = new qx.ui.form.CheckBox( 'use SIFT' );
useSIFTBox.setValue( useSIFT );
useSIFTBox.addListener( "changeValue", e => useSIFT = e.getData() );
container.add( useSIFTBox );

container.add( new qx.ui.basic.Label( 'points/volume:' ) );
const numberOfKeypoints = new qx.ui.form.Spinner( 0, 0, 1000000000 );
numberOfKeypoints.addListener( "changeValue", e =>
    SURF3DParams.numberOfPoints =  RAWParams.numberOfPoints = e.getData() );
numberOfKeypoints.setValue( defaultNumberOfKeypoints );
container.add( numberOfKeypoints );

const useCacheRoot = new qx.ui.form.CheckBox( 'custom cache root' );
container.add( useCacheRoot );
const cacheRoot = new qx.ui.form.TextField("home/temp");
cacheRoot.setVisibility( 'excluded' );
container.add( cacheRoot );
useCacheRoot.addListener( 'changeValue', function ( e ) {

    cacheRoot.setVisibility( e.getData() ? "visible" : "excluded" );

} );

const useCustomPointBox = new qx.ui.form.CheckBox( 'use custom points' );
container.add( useCustomPointBox );
const customPointDirectory = new qx.ui.form.TextField("data/visceral/ruiqi");
customPointDirectory.setVisibility( 'excluded' );
container.add( customPointDirectory );
useCustomPointBox.addListener( 'changeValue', function ( e ) {

    customPointDirectory.setVisibility( e.getData() ? "visible" : "excluded" );

} );


container.add( new qx.ui.basic.Label( "Deformable levels: " ) );
const deformableLevels = new qx.ui.form.Spinner( 0, registrationParams.deformableLevels, 6 );
container.add( deformableLevels );

container.add( new qx.ui.basic.Label( "Iterations/level: " ) );
const deformableIterations = new qx.ui.form.Spinner( 0, registrationParams.deformableIterations, 1000 );
container.add( deformableIterations );

const registerButton = new qx.ui.form.Button( 'register' );
container.add( registerButton );
registerButton.addListener( 'execute', function () {
    computeRegistration().catch( ( e ) => { console.warn( e ); } );
} );

const showBones = new qx.ui.form.CheckBox( 'transform meshes' );
container.add( showBones );

const displayAverageImage = new qx.ui.form.CheckBox( 'Compute average' );
//displayAverageImage.setValue( false );
displayAverageImage.addListener( "changeValue", e =>  {

    for ( let container of [ spacingContainer, displayMeshes, displayTransformedImages ] )
        container.setVisibility( e.getData() ? "visible" : "excluded" );

} );
container.add( displayAverageImage );

const spacingContainer = new qx.ui.container.Composite( new qx.ui.layout.HBox() );
spacingContainer.add( new qx.ui.basic.Label( "spacing (mm)  " ) );
const spacingForm = new qx.ui.form.TextField( spacing.toString() );
spacingForm.addListener( "changeValue", e => spacing = parseFloat( e.getData() ) );
spacingContainer.add( spacingForm, { flex : 1 } );
spacingContainer.setVisibility( "excluded" );
container.add( spacingContainer );

const displayTransformedImages = new qx.ui.form.CheckBox( 'display transformed images' );
container.add( displayTransformedImages );
displayTransformedImages.setVisibility( "excluded" );

const displayMeshes = new qx.ui.form.CheckBox( 'display meshes' );
container.add( displayMeshes );
displayMeshes.setVisibility( "excluded" );

const displayLandmarkBox = new qx.ui.form.CheckBox( 'display landmarks' );
displayLandmarkBox.setValue( displayLandmarks );
displayLandmarkBox.addListener( "changeValue", e => displayLandmarks = e.getData() );
container.add( displayLandmarkBox );

const material = new THREE.MeshBasicMaterial( { color: "red" } );
const material2 = new THREE.MeshBasicMaterial( { color: "black" } );
const geometry = new THREE.PlaneGeometry(  ( 2 + width ) * spacingX, 3 );
geometry.applyMatrix4( new THREE.Matrix4().makeRotationX( 0.5 * Math.PI ) );
/*
for ( let line of [ minLine, maxLine ] ) {

    line.position.x = 0.5 * width * params.spacingX;
//        line.position.y = - 1.5;
    viewer.addMesh( line, { label : "minMaxLine" });

}
*/

const maxLine = new THREE.Line( geometry, material );
const minLine = new THREE.Line( geometry, material );
maxLine.visible = minLine.visible = false;
viewer.addMesh( maxLine );
viewer.addMesh( minLine );

const maxLine2 = new THREE.Line( geometry, material2 );
const minLine2 = new THREE.Line( geometry, material2 );
maxLine2.visible = minLine2.visible = false;
viewer.addMesh( maxLine2 );
viewer.addMesh( minLine2 );


for ( let line of [ minLine, maxLine, minLine2, maxLine2 ] ) {

    line.position.x = 0.5 * width * spacingX;
//        line.position.y = - 1.5;
    viewer.addMesh( line, { label : "minMaxLine" });

}



if ( presets ) {

    if ( presets.organs) organs = presets.organs;

}


async function computeRegistration () {

    const volumesToRegister = meshes.filter( mesh => mesh.visible )
        .map( mesh => ( {

            volume : mesh.userData.file,
            translation : [ 0, 0, 0 ]

        } ) );

    let customPoints;

    if ( useCustomPointBox.getValue() ) {

        const dir = customPointDirectory.getValue();
        customPoints = ( await desk.FileSystem.readDirAsync( dir ) )
            .filter( file => file.name.includes( ".csv" ) )
            .map( file => dir + '/' + file.name );

    }

    const fullMatchParams = Object.assign( {}, matchParams );

    registrationParams.deformableLevels = deformableLevels.getValue();
    registrationParams.deformableIterations = deformableIterations.getValue();

    const reg = new window.FROG.DeformableGroupwiseRegistration( volumesToRegister,
        { SIFT3DParams, SURF3DParams, RAWParams, matchParams : fullMatchParams,
        registrationParams, useRAW, useSIFT, useSURF, customPoints,
        globalParams : useCacheRoot.getValue()? { cacheRoot : cacheRoot.getValue() } : {} } );

    reg.on('log', message => {

        statusLabel.setValue( message );
        console.log( message );

    } );

    let matchStarted = false;
    const numberOfMatches = 0.5 * volumesToRegister.length * ( volumesToRegister.length - 1 );
    let matchCount = 0;
    reg.on('matchLog', message => {

        if ( message.includes( "Pairing..." ) ) {
            matchStarted = true;
            return;
        }

        if ( message.includes( "." ) && matchStarted ) {

            matchCount++;
            statusLabel.setValue( 'Computing pairs : ' +
                Math.round( 100 * matchCount / numberOfMatches ) + '% done ' );

        }

    } );

    let iteration = 1;
    let currentVolume;

    function parseLinear( line ) {
        return line.split( '=' )[ 1 ].split( " " ).map( parseFloat ) ;
    }

    reg.on('registrationLog', message => {

        if ( viewer.isDisposed() ) return;

        for ( let line of message.split( '\n' ) ) {

            if ( line.includes( "translation" ) ) {
                currentVolume = line.match(/\d+/).shift();
                meshes[ currentVolume ].position.fromArray( parseLinear( line ) );
                viewer.render();
            }

            if ( line.includes( "scale" ) ) {
                meshes[ currentVolume ].scale.fromArray( parseLinear( line ) );
                viewer.render();
            }

            if ( line.includes( "average" ) ) {

                statusLabel.setValue( 'Iteration ' + iteration++ + ': ' + line );

            }

        }

    } );

    const output = await reg.execute();
    const registration = output.registration;
    const filefilter2 = getVolumeFilter( filter2 );
    console.log( "registration : ", output );
    console.log(output.volumes);
    output.volumes = output.volumes.filter( vol => filefilter2( vol.volume ) );
    console.log(output.volumes);
    output.spacing = spacing;

    // measure landmark dispersion among volumes from the VISCERAL database
    const stats = await visceral.getLandmarkDistances( output );

    if ( stats.stats ) {

        const glob = stats.global;
        console.log( "stats:", stats.stats );
        console.log( ["mean", "median", "max", "stdev"]
            .map( field => field + ' ' + ( glob[ field ] || 0 ).toFixed( 1 ) )
            .join( ', ' )
        );

        console.log( "all stats:", stats );

    }

    if ( displayLandmarks ) {

        const landmarkViewer = new desk.MeshViewer();
        const radius = 5;
        const geometry = new THREE.SphereGeometry( radius, 32, 32 );
        const scene = landmarkViewer.getScene();

        for ( let name of Object.keys( stats.stats ) ) {

            const rng = new desk.Random( hashCode( name ) );
            const color = [ rng.random(), rng.random(), rng.random() ];
            const material = new THREE.MeshLambertMaterial();
            material.color.fromArray( color );

            for ( let landmark of stats.stats[ name ].positions ) {

                const mesh = new THREE.Mesh( geometry, material );
                mesh.position.fromArray( landmark.position );
                mesh.userData.landmark = landmark;
                landmarkViewer.addMesh( mesh, {  label : name } );

            }

        }

        landmarkViewer.resetView();

        function hashCode ( string ) {

            return string.split("").reduce((a,b) => {a=((a<<5)-a)+b.charCodeAt(0);return a&a;},0);

        }

    }


    let showBonesPromise;
    if ( showBones.getValue() ) showBonesPromise = extractMeshes( output );


    if ( displayAverageImage.getValue() ) await computeAverageImage( output );
    await showBonesPromise;
    statusLabel.setValue( 'All done (' + output.volumes.length + ' volumes)' );

}

viewer.add( new qx.ui.basic.Label( "contrast" ), { bottom : 1, left : 0 } );
const slider = new qx.ui.form.Slider( "horizontal" );
slider.setWidth( 300 );
viewer.add( slider, { left : 60, bottom : 0 } );
slider.setValue( 1 / 0.04 );
slider.addListener( "changeValue", () => {
    viewer.getScene().traverse( object => {
        const slice = object.userData?.viewerProperties?.volumeSlice;
        if ( !slice ) return;
        slice.setContrast( ( slider.getValue() ) * 0.04);
    } );
});



async function computeAverageImage( output ) {

    const commonSpace = new FROG.CommonSpaceMeanImage( output );
    commonSpace.on('log', message => statusLabel.setValue( message ) );
    const result = await commonSpace.execute();
    console.log( result );
    const averageVolume = await result.averageVolume;
    console.log( "Average volume : ", averageVolume );

    if ( !vv ) {

        vv = new desk.VolumeViewer();
        vv.getWindow().addListener( "beforeClose", () => vv = null );

    }

    vv.addVolume( averageVolume, { format : 0, label : filter2 + " (" + output.volumes.length + ' images)' } );

    if ( displayMeshes.getValue() ) {

        const organViewer = new desk.MeshViewer();
        organViewer.getWindow().setCaption( filter1 + " " + filter2 + " (" + output.volumes.length + ' images)' );
        await ACVD.extractMeshes( averageVolume, organs, organViewer );

    }

    if ( displayTransformedImages.getValue() ) {

        const transformedVolumeViewer = new desk.VolumeViewer();

        for ( let [ index, volume ] of result.transformedVolumes.entries() ) {

            await transformedVolumeViewer.addVolume( volume, { label : "" + index } );

        }

    }

}

}

run().catch( ( e ) => {

	console.warn( e );

});


async function extractMeshes( output ) {

    const viewer = new desk.THREE.Viewer();
    viewer.getWindow().setCaption( "All organs" );

    const groups = {};

    for ( let organ of Object.keys( organs ) ) {

        const group = groups[ organ ] = new THREE.Group();
        viewer.addMesh( group, { label : organ } );

    }

    const promises = output.volumes.entries().map( async entry => {

        const [ index, obj ] = entry;
        const meshes = await ACVD.extractMeshes( obj.volume, organs );

        const promises = Object.keys( organs ).map( async organ => {

            const convert = await desk.Actions.executeAsync( {

                action : "mesh2obj",
                inputMesh : meshes[ organ ].file

            } );

            const transform = await desk.Actions.executeAsync( {

                action : "MeshTransform",
                transform : obj.transform,
                source : convert.outputDirectory + "mesh.obj"

            } );

            const mesh = await viewer.addFileAsync( transform.outputDirectory + "output.obj",
                { label : organ + " " + index, parent : groups[ organ ] });

            const opts = organs[ organ ];
			if ( opts.material ) mesh.material.setValues( opts.material );
			if ( opts.transparent != undefined ) mesh.transparent = opts.transparent;
			if ( opts.renderOrder != undefined ) mesh.renderOrder = opts.renderOrder;

        } );

        await Promise.all( promises );

    } );

    await Promise.all( promises );

}


function getVolumeFilter ( filter ) {

    switch ( filter ) {

        case "men":
            return file =>  visceral.getSex( file ) === 1;

        case "women":
            return file =>  visceral.getSex( file ) === 0;

        case "all" :
            return x => true;

        default :
            return file =>  {
                //console.log(visceral.getFileId( file ), filter)
                return visceral.getFileId( file ) === filter;
            };

    }

}
