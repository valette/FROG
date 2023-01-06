"use strict";

const { qx, desk, THREE } = window;

console.clear();

const params = {

    MATCH3D : {
        action : "MATCH3D",
        RansacDist : 20,
        MatchingScale : 1.5,
        MatchingDist : 0.3,
        MatchingDist2 : 0.98,
        ComputeBoundingBoxes : 1,
        RansacMinInliers : 1 //[5, 20, 25, 30, 40]
    },
    SURF3D : {
        action : "SURF3D",
        threshold : 7000,
        target_spacing: 1, 
    },
    MESHING : {
        action : "extract_meshes", 
        threshold : 300,
        max_number_of_vertices : 100000
    },
    inputFiles : [

        "big/00/dea/1.3.12.2.1107.5.1.4.50338.30000012091205384854600003775/1.3.12.2.1107.5.1.4.50338.30000012091205384854600003775.mhd",
        "big/00/dea/1.3.12.2.1107.5.1.4.50338.30000012091205384854600004095/1.3.12.2.1107.5.1.4.50338.30000012091205384854600004095.mhd",
        "big/00/corps//19bad/1.3.46.670589.33.1.32963155271087338575.28270001074118335806/1.3.46.670589.33.1.32963155271087338575.28270001074118335806.mhd"
    ],

    registerAll : false,
    showSlices : true,
    showBoxes : false,
    showRegistration : false
};

const volumes = [];
const meshes = [];
let links;
const boxes = [];
const loadedFiles = [];
const ids = [];
const shuffleButtons = [];
let surfActions = [];
let transform;

const win = new qx.ui.window.Window( "SURF" );
win.set ({layout : new qx.ui.layout.HBox(), width : 700, height : 500});
const tabView = new desk.TabView();
const pane = new qx.ui.splitpane.Pane();
win.add( pane, { flex : 1 } );
if ( desk.auto ) qx.core.Init.getApplication().getRoot().add( pane,
    { width : '100%', height : '100%' } );

function tweakShader (volume) {
    volume.getSlices().forEach(function ( slice, index ) {
        var material = slice.getUserData('mesh').material;
        material.baseShader.extraShaders.push(
            'gl_FragColor[3] *= opacity * step(0.10, rawData[0]);');
        slice.updateMaterial(material);
    } );
}

const filesContainer = new qx.ui.container.Composite();
filesContainer.setLayout( new qx.ui.layout.VBox( 5 ) );
if ( !desk.auto ) tabView.addElement( "files", filesContainer );

const viewer = new desk.SceneContainer( {
    cameraFront : [ 0, -1, 0 ],
    cameraUp : [ 0, 0, -1 ]
});
const statusLabel = new qx.ui.basic.Label( "..." );
viewer.add( statusLabel, { left : "40%", bottom : 10 } );
const tab = tabView.addElement('3D', viewer);
tabView.setSelection( [ tab ] );

const showContainer = new qx.ui.container.Composite();
showContainer.setLayout( new qx.ui.layout.VBox() );
showContainer.setBackgroundColor( "white" );
showContainer.setDecorator( "border-blue" );
viewer.add( showContainer, { bottom : 5, right : 5 } );

const files = [ 0, 1 ].map( ( i ) => {

    const field = new desk.FileField( "" );
    filesContainer.add( new qx.ui.basic.Label( "File " + ( i + 1) ) );
    filesContainer.add( field );
    return field;

} );

const sphereGeometry = new THREE.SphereGeometry( 1, 10, 10 ); 

async function addMesh( file, index ) {

    if ( loadedFiles[ index ] == file ) return;
    if ( meshes[ index ] ) viewer.removeMesh( meshes[ index ], { updateCamera : false } );
    updateInliers();

    async function MPR() {

        viewers[ index ].removeAllVolumes();
        return await viewers[ index ].addVolumeAsync( file );
//        tweakShader( volume );

    }

    const MPRPromise = MPR();

    const meshing = await desk.Actions.executeAsync( {
        input_volume : file,
        ...params.MESHING
    } );

    const mesh = new THREE.Group();
    const updateCamera = index == 0;
    viewer.addMesh( mesh, { label : "Input " + ( index + 1 ), updateCamera } );

    await viewer.addFileAsync( meshing.outputDirectory + "/0.vtk",
        { label : "bones", parent : mesh, updateCamera : index == 0 } );

    meshes[ index ] = mesh;
    setMeshesColor();
    const action = await surfActions[ index ];
    const pointsFile = action.outputDirectory + "pts.json";
    const txt = await desk.FileSystem.readFileAsync( pointsFile, { cache : action } );
    const pts = JSON.parse( txt );
	const matrix = new THREE.Matrix4();
	const material = new THREE.MeshLambertMaterial();
	const imesh = new THREE.InstancedMesh( sphereGeometry, material, pts.points.length );
	const v = new THREE.Vector3();

	for ( let [ i, pt ] of pts.points.entries() ) {

        v.copy( pt );
        matrix.makeScale( pt.scale, pt.scale, pt.scale );
        matrix.setPosition( ...v.toArray() );
		imesh.setMatrixAt( i, matrix );

	}

	imesh.userData.points = pts;
    viewer.addMesh( imesh, { parent : mesh, label : "points", updateCamera } );
    const volume = await MPRPromise;
    viewer.attachVolumeSlices( volume.getSlices(), { parent : mesh, updateCamera } );
    volume.getSlices()[ 0 ].addListener( "changePosition", updatePointsVisibility );
    updatePointsVisibility();
    viewer.render();
    loadedFiles[ index ] = file;
    if ( index == 0 ) viewer.resetView();
    setMeshesColor();

}

async function addMeshesAndMatch() {

	console.log( ids );
	await update( files.map( f => f.getValue() ) );

}

async function match( files ) {

    transform = null;

    surfActions = files.map( file => desk.Actions.executeAsync( {

            input_volume : file,
            outputFileName : "pts",
            writeJSON : 1,
            ...params.SURF3D

        } )

    );

    await Promise.all( surfActions );

    const match = await desk.Actions.executeAsync( {

        input_volume1 : ( await surfActions[ 0 ] ).outputDirectory + "pts.json",
        input_volume2 : ( await surfActions[ 1 ] ).outputDirectory + "pts.json",
        writeInliers : true,
        ...params.MATCH3D

    } );

    const txt = await desk.FileSystem.readFileAsync( match.outputDirectory + "transform.json" );
    const trans = JSON.parse( txt );
    return trans.fail ? null : trans;

}

async function update( files ) {

    try {

        transform = null;
		for ( let button of shuffleButtons ) button.setEnabled( false );
        statusLabel.setValue( "Registering..." );
        const promise = match( files );
        await Promise.all( [ 0, 1 ].map( index => addMesh( files[ index ], index ) ) );
        transform = await promise;
        setMeshesColor();
        updatePointsVisibility();
        statusLabel.setValue( transform ?
			"Done, " + transform.inliers + " inliers." : "Failed");
		for ( let button of shuffleButtons ) button.setEnabled( true );

    } catch( e ) { console.warn( e ); }

}

const layout = new qx.ui.layout.VBox( 5 );
const container = new qx.ui.container.Composite( layout );

const viewers = [ 0, 1 ].map( ( index ) => {

    const MPR = new desk.MPRContainer(null, { nbOrientations : 1 } );
    MPR.maximizeViewer(0);
    container.add(MPR, {flex : 1});
    const fileBrowser = new desk.FileBrowser(null, false);
    fileBrowser.setFileHandler(function () {}); // disable double click
    fileBrowser.setWidth(300);
    fileBrowser.getTree().addListener('changeSelection', function () {

        const selectedFiles = fileBrowser.getSelectedFiles();
        if ( !selectedFiles.length ) return;
        const fileName = selectedFiles[0];
        switch (desk.FileSystem.getFileExtension(fileName)) {
            case "mhd":
            case "png":
            case "jpg":
                files[ index ].setValue(fileName);
                break;
            default :
                break;
        }

    });

    if ( !desk.auto ) tabView.addElement('input ' + ( index + 1 ), fileBrowser );
    return MPR;

} );

pane.add( container );

if ( desk.auto ) pane.add( viewer, 2 );
else {

	pane.add( tabView, 2 );
    win.open();
    win.center();

}

const buttonsContainer = new qx.ui.container.Composite();
buttonsContainer.setLayout( new qx.ui.layout.VBox() );
viewer.add( buttonsContainer, { left : 5, bottom : 5 } );
const switchColors = new qx.ui.form.CheckBox( "switch colors" );
switchColors.addListener("changeValue", setMeshesColor);
showContainer.add(switchColors);

function updatePointsVisibility() {

	for ( let i = 0; i < 2; i++ ) {

		const group = meshes[ i ];
		if ( !group ) continue;
		const arr = group.children.filter( m => m?.userData?.points );
		if (!arr.length ) continue;
		const mesh = arr[ 0 ];
		mesh.instanceMatrix.needsUpdate = true;
		const points = mesh.userData.points;

		let position;
		const arr2 = group.children.filter( m => m?.userData?.viewerProperties?.volumeSlice );

		if ( pointsInSlices.getValue() && arr2.length ) {

			const slice = arr2[ 0 ].userData.viewerProperties.volumeSlice;
			position = slice.getPosition();

		}

        let inliers;
        if ( inlierPoints.getValue() && links ) inliers = links.userData.inliers[ i ];
		const v = new THREE.Vector3();
		const matrix = new THREE.Matrix4();

		for ( let [ i, pt ] of points.points.entries() ) {

			v.copy( pt );
			let scale = pt.scale;

			if ( position != undefined  && ( Math.abs( pt.z - position ) > scale ) )
				scale = 0;

            if ( inliers && !inliers.has( i ) ) scale = 0;
			if ( inlierPoints.getValue() && !transform ) scale = 0;

			matrix.makeScale( scale, scale, scale );
			matrix.setPosition( ...v.toArray() );
			mesh.setMatrixAt( i, matrix );

		}

	}

	viewer.render();

}

function setMeshesColor() {

	pointsInSlices.setVisibility( showPoints.getValue() ? "visible" : "hidden" );
	inlierPoints.setVisibility( showPoints.getValue() ? "visible" : "hidden" );

    for ( let i = 0; i < 2; i++ ) {

        const index = switchColors.getValue() ? 1 - i : i;
        const color = index ? [ 1, 0.7, 0.55 ] : [ 1, 1, 1 ];
        const color2 = index ? [ 1, 0, 0 ] : [ 0, 0, 1 ];
        const group = meshes[ i ];

        if ( group ) {

            const mesh = group.children[ 0 ];

            if ( mesh ) {

                const material = mesh.material;
                material.color.setRGB( ...color );
                material.transparent = true;
                material.opacity = 0.8;
                material.side = THREE.DoubleSide;
                mesh.visible = showMeshes.getValue();
                mesh.renderOrder = 2000;

            }

            const slices = group.children[ 2 ];
            if ( slices ) {
				slices.visible = showSlices.getValue();
				slices.children[ 0 ].renderOrder = -10;
			}
            const points = group.children[ 1 ];

            if ( points ) {
                points.material.color.setRGB( ...color2 );
                points.visible = showPoints.getValue();
            }

        }

        const button = shuffleButtons[ i ];
        if ( button ) button.setDecorator( index ? "button-hover" : "main");

    }

    const mesh = meshes[ 1 ];

    if ( mesh ) {

        if ( transform && !transform.fail && showRegistration.getValue() ) {

            mesh.position.fromArray( transform.translation ).multiplyScalar( -1 );
            mesh.scale.setScalar( 1.0 / transform.scale );

        } else {

            mesh.position.set( 1000, 0, 0 );
            mesh.scale.setScalar( 1.0 );
            if ( ! meshes[ 0 ]?.children[ 0 ] ) return;
            if ( ! meshes[ 1 ]?.children[ 0 ] ) return;
            const box1 = meshes[ 0 ].children[ 0 ].geometry.boundingBox;
            const box2 = meshes[ 1 ].children[ 0 ].geometry.boundingBox;
            const center = new THREE.Vector3();
            box1.getCenter( center );
            box2.getCenter( mesh.position );
            mesh.position.multiplyScalar( -1 ).add( center );
            mesh.position.x += 1000;

        }

    }

    updateInliers();
    viewer.render();

}

function updateInliers() {

    if ( links ) viewer.removeMesh( links );
	updateBoxes();
    links = null;
    if ( !transform ) return;
    const nInliers = transform.inliers;
    const geometry = new THREE.BufferGeometry();
    const material = new THREE.LineBasicMaterial( { color : 0x000000 });
    const positions = new Float32Array( nInliers *  6 );
    const pos = new THREE.BufferAttribute( positions, 3 );
    geometry.setAttribute( 'position', pos );
    const mesh1 = meshes[ 0 ].children[ 1 ];
    const mesh2 = meshes[ 1 ].children[ 1 ];
    const group2 = meshes[ 1 ];
    group2.updateMatrix();
    const matrix1 = new THREE.Matrix4();
    const matrix2 = new THREE.Matrix4();
    const p1 = new THREE.Vector3();
    const p2 = new THREE.Vector3();
    const s = new THREE.Vector3();
    const q = new THREE.Quaternion();
    const inliers = [ new Set(), new Set() ];

    for ( let i = 0; i < nInliers; i ++ ) {

        const [ point1, point2 ] = transform.allInliers[ i ];
        inliers[ 0 ].add( point1 );
        inliers[ 1 ].add( point2 );
        mesh1.getMatrixAt( point1, matrix1 );
        mesh2.getMatrixAt( point2, matrix2 );
        matrix1.decompose( p1, q, s );
        matrix2.decompose( p2, q, s );
        p2.applyMatrix4( group2.matrix );
        pos.setXYZ( i * 2, ...p1.toArray() );
        pos.setXYZ( 1 + i * 2, ...p2.toArray() );

    }
    
    links = new THREE.LineSegments( geometry, material, THREE.LineSegments );
    links.userData.inliers = inliers;
    links.visible = showLinks.getValue();
    viewer.addMesh( links, { label : "links" } );
}

function updateBoxes() {

	for ( let i = 0; i < 2; i++ ) {

		if ( boxes[ i ] ) viewer.removeMesh( boxes[ i ] );
		boxes[ i ] = null;
		if ( !transform || !showBoxes.getValue() ) continue;
		const b = transform[ i ? "bboxB" : "bboxA" ];
		const geometry = new THREE.BoxGeometry();
		const material = new THREE.MeshLambertMaterial( { color : "yellow"} );
		material.opacity = 0.2;
		material.transparent = true;
		material.side = THREE.DoubleSide;
		const box = new THREE.Box3();
		box.min.fromArray( b.min );
		box.max.fromArray( b.max );
		const mesh = new THREE.Mesh( geometry, material );
		mesh.renderOrder = 200000;
		box.getCenter( mesh.position );
		box.getSize( mesh.scale );
		viewer.addMesh( mesh, { label : "box", parent : meshes[ i ] } );
		boxes[ i ] = mesh;

	}

}

const showRegistration = new qx.ui.form.CheckBox('Show registration');
showContainer.add( showRegistration );
showRegistration.setValue( params.showRegistration );
showRegistration.addListener( "changeValue", setMeshesColor );

const showPoints = new qx.ui.form.CheckBox('Show points');
showContainer.add( showPoints );
showPoints.setValue( false );
showPoints.addListener( "changeValue", setMeshesColor );

const pointsInSlices = new qx.ui.form.CheckBox('Points in slices');
showContainer.add( pointsInSlices );
pointsInSlices.setValue( false );
pointsInSlices.setVisibility( "hidden" );
pointsInSlices.addListener( "changeValue", updatePointsVisibility );

const inlierPoints = new qx.ui.form.CheckBox('Inlier points');
showContainer.add( inlierPoints );
inlierPoints.setValue( false );
inlierPoints.setVisibility( "hidden" );
inlierPoints.addListener( "changeValue", updatePointsVisibility );

const showSlices = new qx.ui.form.CheckBox('Show slices');
showContainer.add( showSlices );
showSlices.setValue( params.showSlices );
showSlices.addListener( "changeValue", setMeshesColor );

const showLinks = new qx.ui.form.CheckBox('Show links');
showContainer.add( showLinks );
showLinks.setValue( false );
showLinks.addListener( "changeValue", setMeshesColor );

const showBoxes = new qx.ui.form.CheckBox('Show boxes');
showContainer.add( showBoxes );
showBoxes.setValue( params.showBoxes );
showBoxes.addListener( "changeValue", setMeshesColor );

const showMeshes = new qx.ui.form.CheckBox('Show meshes');
showContainer.add( showMeshes );
showMeshes.setValue( true );
showMeshes.addListener( "changeValue", setMeshesColor );

let windowIsClosed;

win.addListener('close', function () {

    windowIsClosed = true;
    for ( let MPR of viewers ) MPR.dispose();
    tabView.dispose();
    win.destroy();
    viewer.dispose();

} );

desk.FileSystem.traverse( "big/00/corps/", traverse, afterTraverse );

function traverse( file ) {

    if ( desk.FileSystem.getFileExtension(file) === "mhd" ) volumes.push (file);

}

async function afterTraverse() {

    for ( let file of params.inputFiles.slice( 0,2 ) )
        if ( !volumes.includes( file ) ) volumes.push( file );

    console.log( volumes.length + " files" );
    volumes.sort( ( a, b ) => a.localeCompare( b ) );
    addRandomButton( 0 );
    addRandomButton( 1 );
    setMeshesColor();


    if ( desk.auto || !params.registerAll ){

		for ( let [ i, field ] of files.entries() ) {
			field.setValue( params.inputFiles[ i ] );
			ids[ i ] = volumes.indexOf( params.inputFiles[ i ] );
		}

        addMeshesAndMatch();

		for ( let field of files )
			field.addListener( "changeValue", addMeshesAndMatch );

		return;

    }

    for ( let volume of volumes ) {

        for ( let volume2 of volumes ) {

            if ( windowIsClosed ) break;
            await update( [ volume, volume2 ] );

        }

    }

    console.log( "done" );

}

viewer.add( new qx.ui.basic.Label( "contrast" ), { top : 10, right : 0 } );
const slider = new qx.ui.form.Slider( "vertical" );
slider.setHeight( 300 );
viewer.add( slider, { right : 10, top : 30 } );
slider.setValue( 100 - 1 / 0.04 );
slider.addListener( "changeValue", () => {
    viewer.getScene().traverse( object => {
        const slice = object.userData?.viewerProperties?.volumeSlice;
        if ( !slice ) return;
        slice.setContrast( ( 100 - slider.getValue() ) * 0.04);
    } );
});

function addRandomButton( index ) {

	const button = new qx.ui.form.Button( "New " + ( index + 1 ) );
	const rng = new desk.Random( index * 1000 );
	shuffleButtons[ index ] = button;

	button.addListener("execute", function () {

		let id = ids[ index ];
		while ( id == ids[ index ] ) id = Math.floor( rng.random() * volumes.length );
		ids[ index ] = id;
		files[ index ].setValue( volumes[ id ] );

    } );

    buttonsContainer.add( button );

}

