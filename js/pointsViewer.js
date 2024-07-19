'use strict';
console.clear();
const { async, THREE, desk, qx } = window;

const params = {
    volume : "big/visceral/volumes/10000100_1_CTce_ThAb.nii.gz",
    spacing : 0.75,
    numberOfPoints : 20000,
    visibilityRadius : 5
};

const viewer = new desk.THREE.Viewer();
const MPR = new desk.MPR.Container();
viewer.getWindow().add( MPR, { flex : 1 } );

let volume;
let volumeMeshes;
let pointsMesh;
let points;

const label = new qx.ui.basic.Label( "" );
viewer.add( label, { right : 10, bottom : 10 } );

const surf = new desk.Action( "SURF3D");
surf.setParameters( { type : 3, writeJSON : 1, writeCSVGZ : 0, writeCSV : 0, writeBIN : 0, outputFileName : "points"}, true );
surf.setParameters( { spacing : params.spacing, numberOfPoints : params.numberOfPoints } );
surf.getForm( "inputVolume").addListener( "changeValue", updateVolume );
surf.addListener( "actionUpdated", updatePoints );
viewer.getWindow().addAt( surf, 0 );

const viewerVisible = new qx.ui.form.ToggleButton( "hide slices" );
viewerVisible.addListener( "changeValue", e => MPR.setVisibility( !e.getData() ? "visible" : "excluded" ) );
surf.add( viewerVisible );

async function updateVolume() {

    if ( volumeMeshes ) viewer.removeMeshes( volumeMeshes );
    volumeMeshes = null;
    if ( pointsMesh ) viewer.removeMeshes( pointsMesh );
    pointsMesh = null;
    MPR.removeAllVolumes();
    const newFile = surf.getForm( "inputVolume").getValue();
    volume = await MPR.addVolumeAsync( newFile );
    const slices = volume.getSlices();
    for ( let slice of slices ) slice.addListener( "changePosition", updatePointsVisibility );
    volumeMeshes = viewer.attachVolumeSlices( slices );
    
}


async function updatePoints( event ) {

    if ( pointsMesh ) viewer.removeMesh( pointsMesh );
    pointsMesh = null;

    const txt = await desk.FileSystem.readFileAsync( surf.getOutputDirectory() + "points.json",
        { cache : event.getData().response } );

    points = JSON.parse( txt ).points;
    const geometry = new THREE.SphereGeometry( 1, 32, 16 ); 
    const material = new THREE.MeshLambertMaterial( { color: 0xffffff } ); 
    pointsMesh = new THREE.InstancedMesh( geometry, material, points.length );
    viewer.addMesh( pointsMesh );
    label.setValue( points.length + " points" );
    updatePointsVisibility();

}

function updatePointsVisibility() {

    if ( !pointsMesh ) return;
    const slices = volume.getSlices();
    const pos = [ 1, 2, 0 ].map( i => slices[ i ].getPosition() );
    const matrix = new THREE.Matrix4();
    const colors = [ 0xCD1719, 0xFFED00, 0x009FE3  ].map( c => new THREE.Color( c ) );
    let color = colors[ 0 ];

	for ( let i = 0; i < points.length; i ++ ) {

        const pt = points[ i ];
        const xyz = [ pt.x, pt.y, pt.z ];
        let visible = false;
        for ( let j = 0; j < 3; j++)
            if ( Math.abs( xyz[ j ] - pos[ j ] ) < params.visibilityRadius ) {
                visible = true;
                color = colors[ j ];
            }
        
        const scale = visible ? pt.scale : 0;
        matrix.makeScale( scale, scale, scale );
        matrix.setPosition(  points[ i ].x, points[ i ].y, points[ i ].z );
		pointsMesh.setColorAt( i, color );
		pointsMesh.setMatrixAt( i, matrix );
	}

    pointsMesh.instanceMatrix.needsUpdate = true;
    pointsMesh.instanceColor.needsUpdate = true;
    viewer.render();
}

surf.setParameters( { inputVolume : params.volume } );
