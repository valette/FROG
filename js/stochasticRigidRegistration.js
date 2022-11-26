'use strict';

const params = {

    maxNumberOfIterations : 200,
    batchSize : 50, // stochastic batch size
    maxSize : 10000, // maximum number of volumes
    learningRate : 0.2, // learning rate
    sortStart : 3, // time at which volumes are displayed in z-sorted order
    outlierStart : 5, // time at which removing outlier starts
    outlierRatio : 10, // outlier threshold : displacement > outlierRatio * medianDisplacement
    spacingX : 500, // grid x spacing for display
    spacingZ : 3000, // grid z spacing for display
    animationDuration : 500, // length of animation in milliseconds
    showDisplacementNorms : false,
    displayMeanZ : false,
    loadVolumeConcurrency : 5

};

const { qx, desk, THREE, _, async, alert } = window;
desk.URLParameters.parseParameters( params );
console.clear();

( async function () {

const placesAll = [
    'big/01_Registrations/registration-001-931.json',
    //'big/01_MHD',
    //'big/01_MHD',
    //'big/00_MHD',
    'big/visceral/volumes',
    'big/visceral/silver/volumes/CTce_ThAb',
    'big/visceral/volumes_CT_wb'
];

const places = placesAll;
const outliers = []; // this array will contain detected outliers
let files = []; // this array will be filled with volume files to register
const lines = {};
let width;

// setup UI
const viewer = new desk.MeshViewer( {

    orthographic : true,
    cameraFront : [ 0, 1, 0 ],
    cameraUp : [ 0, 0, 1 ]

} );

if ( desk.auto ) viewer.getOptionsButton().setVisibility( "excluded" );
if ( desk.auto ) viewer.fillScreen();

const resetView = _.throttle( () => viewer.resetView(), 500 );
viewer.getControls().noRotate = true;
const statusLabel = new qx.ui.basic.Label( '' );
statusLabel.set( { backgroundColor : "white" } );
statusLabel.setValue = _.throttle( statusLabel.setValue.bind( statusLabel ), 500 );
viewer.add( statusLabel, { left : "30%", bottom : 50 } );
const rng = new desk.Random( "random" );
const batchSizeLabel = new qx.ui.form.Spinner( 2, 2, 300 );
viewer.add( batchSizeLabel, { left : "81%", bottom : 10 } );
viewer.add( new qx.ui.basic.Label(' batch size '), { right : "20%", bottom : 10 } );

batchSizeLabel.addListener( 'changeValue', function ( ) {

    width = Math.max ( 4 * Math.round( Math.sqrt( size ) ), batchSizeLabel.getValue() );

} );


function addFile ( file ) {

    if ( !file.endsWith( ".mhd" ) && !file.endsWith( ".nii.gz" ) ) return;
    files.push( file );

}

// build file list
for ( let place of places ) {

    if ( place.endsWith( ".json" ) ) {

        const content = JSON.parse( await desk.FileSystem.readFileAsync( place ) );
        const list = Array.isArray( content ) ? content : content.volumes;
        files.push( ...list );

    } else
        await desk.FileSystem.traverseAsync( place, addFile );

}

files = _.uniq( files );
files.length = Math.min( files.length, params.maxSize );
const size = files.length;
batchSizeLabel.setValue( params.batchSize );

const meshes = await async.mapLimit( files, params.loadVolumeConcurrency, async file => {

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
    const slice = mesh.children[ 0 ];
    slice.renderOrder = 2000;
    slice.material.depthTest = false;
    slice.material.transparent = false;
    anchor.position.x = ( id % width ) * params.spacingX;
    anchor.position.z = - Math.floor( id / width ) * params.spacingZ;
    anchor.userData.initial = anchor.position.clone();
    anchor.userData.final = new THREE.Vector3();
    statusLabel.setValue( ( id + 1 ) + " volumes loaded" );
    resetView();
    return mesh;

} );

viewer.resetView();
statusLabel.setValue( "All " + meshes.length + " volumes loaded. Ready to iterate" );
const button = new qx.ui.form.Button( 'iterate' );
const spinner = new qx.ui.form.Spinner( 0, 1, 1e10 );

viewer.add( button, { bottom : 10, left : "50%" } );
viewer.add( spinner, { bottom : 10, right : "50%" } );

let numberOfIterations = 0;

button.addListener( 'execute', async function () {

    try {

        button.setEnabled( false );
        const defaultSpinnerValue = spinner.getValue();
    
        while ( spinner.getValue() ) {
    
            await iterate();
            spinner.setValue( spinner.getValue() - 1 );
    
        }
    
        button.setEnabled( true );
        spinner.setValue( defaultSpinnerValue );
        if ( numberOfIterations >= params.maxNumberOfIterations )
            for ( let o of [ button, spinner ] ) o.setVisibility( "excluded" );

    } catch( e ) { console.warn( e ); }

} );

if ( !desk.auto ) addFilterButton();

function addFilterButton() {

    const button = new qx.ui.form.ToggleButton( 'filter' );
    viewer.add( button, { top : 10, right : 10 } );

    const container = new qx.ui.container.Composite( new qx.ui.layout.VBox() );
    container.setVisibility( 'excluded' );
    container.setWidth( 100 );
    button.addListener( 'changeValue', function ( e ) {

        if ( e.getData() ) {

            container.setVisibility( 'visible' );
            filter();
            maxLine.visible = minLine.visible = true;

        } else {

            container.setVisibility( 'excluded' );
            meshes.forEach( mesh => { mesh.visible = true; } );
            maxLine.visible = minLine.visible = false;

        }

        viewer.render();

    } );

    container.addListenerOnce( 'appear', function () {

        const z = _.mean( meshes.map( mesh => mesh.userData.z ) );
        if ( isNaN( z ) ) return;
        max.setValue( z );
        min.setValue( z );

    } );

    viewer.add( container, { top : 40, right : 10 } );

    const max = new qx.ui.form.Spinner( -10000, 0, 10000 );
    container.add( max );

    const min = new qx.ui.form.Spinner( -10000, 0, 10000 );
    container.add( min );

    const label = new qx.ui.basic.Label( '' );
    container.add( label );

    const exportButton = new qx.ui.form.Button( 'export' );
    container.add( exportButton );
    exportButton.addListener( 'execute', saveRegistration );

    min.addListener( 'changeValue', filter );
    max.addListener( 'changeValue', filter );

    const material = new THREE.MeshBasicMaterial( { color: "red" } );
    const geometry = new THREE.PlaneGeometry(  ( 2 + width ) * params.spacingX, 3 );
    geometry.applyMatrix4( new THREE.Matrix4().makeRotationX( 0.5 * Math.PI ) );
    const maxLine = new THREE.Mesh( geometry, material );
    const minLine = new THREE.Mesh( geometry, material );

    for ( let line of [ minLine, maxLine ] ) {

        line.position.x = 0.5 * width * params.spacingX;
        line.position.y = - 1.5;
        viewer.addMesh( line, { label : "minMaxLine" });

    }

    maxLine.visible = minLine.visible = false;

    function filter() {

        let n = 0;
        const zMin = min.getValue();
        const zMax = max.getValue();

        minLine.position.z = zMin;
        maxLine.position.z = zMax;

        for ( let mesh of meshes ) {


            if ( !mesh.userData.outlier &&
                ( ( mesh.userData.zMax + mesh.position.z ) >= zMax ) &&
                ( ( mesh.userData.zMin + mesh.position.z ) <= zMin ) ) {

                mesh.visible = true;
                const i = n % width;
                const j = Math.floor( n / width );
                mesh.parent.position.set( i * params.spacingX, 0, - j * params.spacingZ );
                n++;

            } else {

                mesh.visible = false;

            }

        }

        label.setValue( n + ' volumes' );
        viewer.render();

    }

    async function saveRegistration () {

        var obj = {
            volumes : [],
            positions : [],
            zMin : min.getValue(),
            zMax : max.getValue()
        };

        for ( let mesh of meshes.filter( mesh => mesh.visible ) ) {

                obj.volumes.push( mesh.userData.file );
                obj.positions.push( 
                    mesh.position.x,
                    mesh.position.y,
                    mesh.position.z,
                    0
                );

        }

        const file = 'data/stochasticRegistration.json';
        await desk.FileSystem.writeJSONAsync( file, obj );
        alert ( "file written : " + file );

    }

}

async function iterate () {

    numberOfIterations++;

    const batch = getRandomVolumes();

    setupAnimation( batch );

    // construct file list on disk
    const files = batch.map( mesh => mesh.userData.viewerProperties.label );
    const fileOnDisk = await desk.FileSystem.writeCachedFileAsync ( "files.json",
        JSON.stringify( files ) );

    // launch registration
    const registrationDone = desk.Actions.executeAsync( {

        action : "LSRegistration",
        files : fileOnDisk

    } );

    await animate();

    // read registration results
    const registration = JSON.parse( await desk.FileSystem.readFileAsync(
            ( await registrationDone ).outputDirectory + "registration.json") );

    const currentMean = batch.reduce( ( prev, curr ) => prev.add( curr.position ),
        new THREE.Vector3() ).multiplyScalar( 1 / batch.length );

    const newPositions = batch.map( ( mesh, index ) => {

        return new THREE.Vector3().fromArray( registration.positions, 4 * index );

    });

    const newMean = newPositions.reduce( ( prev, curr ) => prev.add( curr ),
        new THREE.Vector3() ).multiplyScalar( 1 / batch.length );

    batch.forEach( ( mesh, index ) => {

        const alpha = mesh.userData.registeredOnce ? params.learningRate : 1;

        const initial = mesh.position.clone().sub( currentMean );
        const final = newPositions[ index ].clone().sub( newMean );
        const displacement = initial.clone().lerp( final, alpha ).sub( initial );

        mesh.userData.registeredOnce = true;

        // displace mesh according to registration
        mesh.position.add( displacement );
        mesh.userData.displacement = displacement.length() / alpha;

        // displace slice frame to give a hint on convergence status
        const frame = mesh.children[ 0 ].children[ 0 ];
        frame.position.copy( displacement ).multiplyScalar( 1 / alpha );
        const material = frame.material;
        material.color.set( 'blue' );
        material.opacity = 0.5;
        material.transparent = false;

        if ( params.showDisplacementNorms ) {

            // display displacement norm on each mesh
            const sprite = mesh.userData.sprite;

            if ( sprite ) {

                mesh.parent.remove( sprite );
                sprite.material.dispose();

            }

            mesh.userData.sprite = new THREE.TextSpriteHelper( mesh.parent,
                '' + Math.round( displacement.length() / alpha ) );

        }

    });

    computeStatistics();
    viewer.render();
    
}

function computeStatistics() {

    const registeredMeshes = meshes.filter ( mesh => ( !mesh.userData.outlier && mesh.userData.registeredOnce ) );
    if ( !registeredMeshes.length ) return;

    const displacements = registeredMeshes
        .map ( mesh => mesh.userData.displacement )
        .sort();

    const mean = _.mean( displacements );
    const median = displacements[ Math.round( displacements.length / 2 ) ];
    const maxMesh = _.maxBy( registeredMeshes, mesh => mesh.userData.displacement );
    const max = maxMesh.userData.displacement;

    for ( let mesh of registeredMeshes )
        mesh.children[ 0 ].children[ 0 ].material.color.set( 'blue' );

    maxMesh.children[ 0 ].children[ 0 ].material.color.set( 'red' );
    const batchSize = params.batchSize;

    statusLabel.setValue(
        meshes.length + " volumes, " +
        "Iteration " + numberOfIterations +
        ' ( ' + ( ( numberOfIterations * batchSize ) / size ).toFixed(1) +
        ' / vol, gain = ' +
        ( size * size / ( numberOfIterations * batchSize * batchSize ) ).toFixed( 1 ) + ' ), ' +
        "Mean = " + mean.toFixed( 1 ) + ", " +
        " Median = " + median.toFixed( 1 ) + ", " +
        " Max = " + max.toFixed( 1 ) 
    );

    // possibly remove one outlier
    if  ( ( numberOfIterations >= params.outlierStart * size / batchSizeLabel.getValue() ) &&
        ( max > ( params.outlierRatio * median ) ) ) {

        outliers.push( maxMesh );
        maxMesh.userData.outlier = true;
        shuffledArray.length = 0;

    }

}

async function animate() {

    const start = performance.now();

    while ( true ) {

        let ratio = Math.min( 1, ( performance.now() - start ) / params.animationDuration );

        for ( let mesh of meshes ) {

            const anchor = mesh.parent;
            const pos = anchor.userData;
            anchor.position.lerpVectors( pos.initial, pos.final, ratio ); 

        }

        await viewer.renderAsync();
        if ( ratio >= 1 ) break;

    }

    for ( let mesh of meshes ) {

        const pos = mesh.parent.userData;
        pos.initial.copy( pos.final );

    }

}

function setupAnimation ( batch ) {

    if ( numberOfIterations < params.sortStart * size / batchSizeLabel.getValue() ) {

        // setup animation to slide all volumes down
        for ( let mesh of meshes ) {
    
            const anchor = mesh.parent;
            anchor.userData.final.copy( anchor.position );
            anchor.userData.final.z = anchor.position.z - params.spacingZ;
    
        }

    } else sortMeshes();

    // setup animation for next volume batch
    batch.forEach( ( mesh, index ) => {

        const finalPos = mesh.parent.userData.final;
        finalPos.z = params.spacingZ;
        finalPos.x = index * params.spacingX;

    } );

    // put rejected meshes aside
    outliers.forEach( ( mesh, index ) => {

        const pos = mesh.parent.userData;
        pos.final.z = - params.spacingZ * Math.round( index / 5 );
        pos.final.x = - params.spacingX * ( ( index % 5 ) + 4 );

    } );

}

function sortMeshes() {

    // setup animation to sort meshes;
    meshes.forEach( mesh => {

        const bounds = mesh.children[ 0 ].userData.viewerProperties.volumeSlice.getBounds();
        mesh.userData.zMin = bounds [ 4 ];
        mesh.userData.zMax = bounds [ 5 ];
        const z = mesh.userData.z = 0.5 * ( bounds [ 4 ] + bounds [ 5 ] ) + mesh.position.z;

        if ( params.displayMeanZ ) {

            // display mean z on each mesh
            var sprite = mesh.userData.Zsprite;

            if ( sprite ) {

                mesh.parent.remove( sprite );
                sprite.material.dispose();

            }

            mesh.userData.Zsprite = new THREE.TextSpriteHelper( mesh.parent, '' + Math.round( z ) );

        }

    } );

    const sortedMeshes = _.sortBy( meshes.filter( mesh => !mesh.userData.outlier && mesh.userData.registeredOnce ),
        mesh => - mesh.userData.z );

    sortedMeshes.forEach( ( mesh, index ) => {

        const finalPos = mesh.parent.userData.final;
        const i = index % width;
        const j = Math.floor( index / width );
        finalPos.x =  i * params.spacingX;
        finalPos.z = - params.spacingZ * ( 1 + j );

        let line = lines[ j ];

        if ( !line ) {

            const material = new THREE.MeshBasicMaterial( { color: 0x00dddd} );
            const geometry = new THREE.PlaneGeometry(  ( 2 + width ) * params.spacingX, 10 );
            geometry.applyMatrix4( new THREE.Matrix4().makeRotationX( 0.5 * Math.PI ) );
            geometry.computeBoundingBox();
            geometry.computeBoundingSphere();
            line = lines[ j ] = new THREE.Mesh( geometry, material );
            line.position.x = 0.5 * width * params.spacingX;
            viewer.getScene().add( line );
            line.userData.iteration = -1;

        }

        const data = line.userData;

        if ( data.iteration != numberOfIterations ) {

            data.iteration = numberOfIterations;
            data.nZ = 0;
            data.sZ = 0;

        }

        data.nZ++;
        data.sZ += mesh.userData.z;

    } );

    for ( let j of Object.keys( lines ).map( parseFloat ) ) {

        const line = lines[ j ];
        const data = lines[ 0 ].userData;
        line.position.z = - ( j + 1 ) * params.spacingZ + data.sZ / data.nZ;

    }

}

let shuffledArray;

function getRandomVolumes() {

    const vols = [];

    while ( vols.length < batchSizeLabel.getValue() ) {

        if ( !shuffledArray || !shuffledArray.length ) {

            shuffledArray = meshes.slice();
            const l = shuffledArray.length;

            for ( let j = 0; j < l * 4; j++ ) {

                const id1 = Math.round( ( l - 1 ) * rng.random() );
                const id2 = Math.round( ( l - 1 ) * rng.random() );
                const temp = shuffledArray[ id1 ];
                shuffledArray[ id1 ] = shuffledArray[ id2 ];
                shuffledArray[ id2 ] = temp;

            }

        }

        const mesh = shuffledArray.pop();
        if ( !mesh.userData.outlier ) vols.push( mesh );

    }

    return vols;

}

viewer.add( new qx.ui.basic.Label( "contrast" ), { top : 80, right : 0 } );
const slider = new qx.ui.form.Slider( "vertical" );
slider.setHeight( 300 );
viewer.add( slider, { right : 10, top : 100 } );
slider.setValue( 100 - 1 / 0.04 );
slider.addListener( "changeValue", () => {
    viewer.getScene().traverse( object => {
        const slice = object.userData?.viewerProperties?.volumeSlice;
        if ( !slice ) return;
        slice.setContrast( ( 100 - slider.getValue() ) * 0.04);
    } );
});

} )  ();
