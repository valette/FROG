"use strict";

const { require, MHD, alert, async, desk, qx, THREE, _, operative, numeric, laplaceSolver } = window;

console.clear();
const solverWorker = laplaceSolver.worker();
var dataDir = "big/00";
//var dataDir = "big/01_MHD";
//dataDir = "big/visceral/volumes";

var defaultFile = 'data/volumes-01_MHD-200.json';
//defaultFile = 'data/visceral/10000109_29193.json';
//defaultFile = 'data/visceral/10000109_187.json';
defaultFile = 'data/volumes-01_MHD-10.json';

const loadDefaultFile = true;

const opts = {
    orthographic : false,
    spacing : 500,
    edgeRemovalRatio : 0.04,
    finalEdgesRatio : 7,
    useScale : true,
    contrast : 2,
    displayInliers : true,
    hideImageBackgrounds : false
};

const params = {
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
        target_spacing: 1.5,
        writeJSON : 1,
        outputFileName : "pts"
//        maximum_dimension :100
    }

};

desk.URLParameters.parseParameters( opts );
var version = "17";
var useCache = true;
var dataSetSize = 20;
var cleanDB = true;
var snapshots = false;
var saveCache = true;
var textureSamplingRatio;

var files;

var width = Math.round(Math.sqrt(dataSetSize));

String.prototype.hashCode = function() {
    var i, l, hval =  0x811c9dc5;

    for (i = 0, l = this.length; i < l; i++) {
        hval ^= this.charCodeAt(i);
        hval += (hval << 1) + (hval << 4) + (hval << 7) + (hval << 8) + (hval << 24);
    }
    return ("0000000" + (hval >>> 0).toString(16)).substr(-8);
};


function getHash() {
    return (JSON.stringify(volumes) + version + JSON.stringify(params)).hashCode();
}

var dummyVector = new THREE.Vector3();
function getPosition(index) {
    return dummyVector.set((index % width) * opts.spacing, Math.floor(index / width) * opts.spacing, 0);
}


var volumes = [];
var matches = [];
var pointsFile = [];


var positions = [];
var meshes = [];

var valence;

var timings = [];
var startTime;

async function surf3d( volume ) {

    const surf3D = await desk.Actions.executeAsync(
        Object.assign( {
            action : "SURF3D", 
            input_volume:volume
        }, params.surf3d ) );
    
    const point3dfile = surf3D.outputDirectory + "pts.json";

    await Promise.all( pointsFile.map( async function( filename ) {
        const match = await desk.Actions.executeAsync(
            _.merge({action : "MATCH3D", 
                input_volume1 : filename, 
                input_volume2 : point3dfile,
                stdout: true,
//                silent : true
        }, params.match3d) );

        matches[ pointsFile.indexOf( filename ) ][ pointsFile.length ] = JSON.parse( match.stderr );
    } ) );

    timings.push((performance.now() - startTime) / 1000);
    updateStatus("Volume " + pointsFile.length + "/" + volumes.length + " done");
    pointsFile.push(point3dfile);
}


function getMatrice(i, j){
    return matches[Math.min(i, j)][Math.max(i, j)];
}

function resetPositions() {

    positions = volumes.map(() => [0, 0, 0, 0] );
    displayAll();

}


async function startRegistration ( ) {

    const beginTime = performance.now();
    let iteration = 0;

    while ( 1 ) {

        const state = await solverWorker.iterate( true );
        let ok;

        if (state.error) {
            ok = false;
            alert("error : " + state.error);
            return;
        }

        const numberOfEdges = state.numberOfEdges;
        valence = state.valences;
        updatePositions(state.positions);
        ok = numberOfEdges > matches.length * opts.finalEdgesRatio;
        const message = "after " + Math.round((performance.now() - beginTime) / 1000) +
            "s. : " + ( iteration + 1 )+ " iterations : " + numberOfEdges + " edges" +
            (!ok ? "  ***finished***" : "");

        if ( snapshots ) {
            viewer.snapshot( {
                path : "data/snapshots/snap" + ( iteration + 1000 ) + ".png",
                ratio : 3
            } );
        }

        iteration++;
        displayStatus( message );

        if ( !removeEdgesButton.getValue() || !ok ) {
            removeEdgesButton.setValue(false);
            return;
        }


    }
}

function updatePositions( volumesPositions ) {

    const z3 = 0;
    for (let i = 0 ; i < matches.length ; i++) {
        let pos = positions[i] = positions[i] || [];
        let i4 = 4 * i;
        pos[0] = volumesPositions[i4] - volumesPositions[z3];
        pos[1] = volumesPositions[i4 + 1] - volumesPositions[z3 + 1];
        pos[2] = volumesPositions[i4 + 2] - volumesPositions[z3 + 2];
        pos[3] = volumesPositions[i4 + 3] - volumesPositions[z3 + 3];
    }

    displayAll();

}

async function loadAll( ){

//viewer.setViewpoint(JSON.parse('{"controls":{"zoomStart":0.15897047691143074,"zoomEnd":0.15897047691143074,"panStart":{"x":0.48485995457986375,"y":0.22558667676003027},"panEnd":{"x":0.48485995457986375,"y":0.22558667676003027},"eye":{"x":90.94013971634195,"y":-3262.4301040152654,"z":436.929711843971},"objectUp":{"x":0.0025487766324721658,"y":0.13281186357592506,"z":0.9911379886933694},"objectPosition":{"x":3229.878654210379,"y":-3214.737577962988,"z":30.180340145935418},"target":{"x":3138.938514494037,"y":47.6925260522775,"z":-406.7493716980356}},"camera":{"metadata":{"version":4.3,"type":"Object","generator":"ObjectExporter"},"object":{"uuid":"94C80EDA-86F2-4C0C-9C64-DE1069DF7EDC","type":"PerspectiveCamera","fov":60,"aspect":1.7463617463617465,"near":3.360015168602809,"far":33600.15168602809,"matrix":[0.9996153116226196,0.027034802362322807,-0.006193222478032112,0,0.0025487763341516256,0.13281187415122986,0.9911379814147949,0,0.02761775255203247,-0.9907724857330322,0.1326918751001358,0,3229.878662109375,-3214.737548828125,30.180339813232422,1],"children":[{"uuid":"9068CC52-93EB-4CB2-8A8C-F1BFD41A3AAC","type":"DirectionalLight","color":11184810,"intensity":1,"matrix":[1,0,0,0,0,1,0,0,0,0,1,0,0.19245009124279022,0.19245009124279022,0.9622504711151123,1]},{"uuid":"D67E3D99-8949-42D0-8292-1D7661F9AC13","type":"Object3D","matrix":[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]}]}},"bbdl":3360.015168602809}'));
    const v2 = {"controls":{"zoomStart":0.10868245294474803,"zoomEnd":0.10868245294474803,"panStart":{"x":0.4286581663630844,"y":0.12264723740133576},"panEnd":{"x":0.4286581663630844,"y":0.12264723740133576},"eye":{"x":-40.48208860847251,"y":-733.0352602748914,"z":69.39274413311932},"objectUp":{"x":0.01398132146440316,"y":0.0934691048776026,"z":0.9955240072863395},"objectPosition":{"x":723.5508443076378,"y":-758.1060423183898,"z":-386.3598443681827},"target":{"x":764.0329329161103,"y":-25.070782043498426,"z":-455.752588501302}},"camera":{"object":{"uuid":"A706103E-5FB0-4414-912D-3CD45E8A8F5B","type":"PerspectiveCamera","matrix":[0.9983941316604614,-0.055966537445783615,-0.008766969665884972,0,0.013981322757899761,0.09346909821033478,0.995523989200592,0,-0.05489658936858177,-0.9940479397773743,0.09410148859024048,0,723.5508422851562,-758.1060180664062,-386.3598327636719,1],"children":[{"uuid":"999988B8-85B2-4593-A69D-23D530F3C574","type":"DirectionalLight","matrix":[1,0,0,0,0,1,0,0,0,0,1,0,0.19245009124279022,0.19245009124279022,0.9622504711151123,1],"color":11184810,"intensity":1},{"uuid":"3F6D7D30-F405-43B7-8E6D-E9389A4202AE","type":"Object3D","matrix":[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]}],"zoom":1,"fov":60,"aspect":3.21227621483376,"near":1.1161121217207792,"far":11161.121217207792}},"bbdl":1116.1121217207792};
//	viewer.setViewpoint(v2);
	let count = 0;

    await async.eachLimit( volumes, 8, async function ( volume ) {

        const index = volumes.indexOf( volume );

        const mesh = await viewer.addVolumeAsync(
            volume, { sliceWith : { sampling_ratio : textureSamplingRatio }, position : getPosition( index ) }
        );

        for ( let slice of mesh.children ) {
            if ( opts.hideImageBackgrounds ) {
                slice.material.baseShader.extraShaders.push('if (rawData[0] < 0.10) discard;');
                slice.children[ 0 ].visible = false;
                slice.userData.viewerProperties.volumeSlice.updateMaterial(slice.material);
            }
            slice.material.transparent = false;
            slice.userData.viewerProperties.volumeSlice.setBrightnessAndContrast( 0, opts.contrast );
        }

        count++;
		displayStatus("Loaded " + count + "/" + volumes.length + " volumes");
        meshes[index] = mesh;
        resetView();
     } );
}

function displayAll(){

    for ( let [ id, mesh ] of meshes.entries() ) {

        const translation = positions[ id ];
        const scale = opts.useScale ? Math.exp( translation[ 3 ] ) : 1;
        mesh.scale.setScalar( scale );
        mesh.position.fromArray(translation);
        mesh.position.add(getPosition(id));

    }

    updateSprites();
    viewer.render();

}

function updateSprites(){
    if (!opts.displayInliers || typeof valence === "undefined") {
        return;
    }

    const bbox = new THREE.Box3();

    meshes.forEach( function ( mesh, id ) {
        var texture, canvas, context;
        var sprite = mesh.userData.sprite;
        var height = 32;
        if (!sprite) {
            canvas = document.createElement('canvas');
            canvas.height = height;

            context = canvas.getContext("2d");
            var width = context.measureText(volumes.length + '').width * height / 10;
            canvas.width = Math.pow(2, Math.ceil(Math.log(width) / Math.log(2)));
            texture = new THREE.Texture(canvas);
            var material = new THREE.SpriteMaterial({map :texture});//, depthTest : false});
            mesh.userData.sprite = sprite = new THREE.Sprite(material);
        
            bbox.setFromObject(mesh);
            bbox.min.sub(mesh.position);
            bbox.max.sub(mesh.position);
            sprite.position.copy(bbox.min).add(bbox.max).multiplyScalar(0.5);
            sprite.position.z = bbox.max.z;
            var size = 100;
            sprite.scale.x = size * width/ height;
            sprite.scale.y = size;
            mesh.add(sprite);
        }
        sprite.visible = showSprites.getValue();
        texture = sprite.material.map;
        canvas = texture.image;
        context = canvas.getContext("2d");
        var text = valence[id] + '';
        context.clearRect(0, 0, canvas.width, canvas.height);
        context.fillStyle = 'rgba(0, 0, 0, 0.5)';
        context.fillRect(0, 0, canvas.width, canvas.height);
//        context.fillText(text, 0, height - 5);
        context.font = height + 'px bold Arial';
        context.fillStyle = "yellow";
        context.fillText(text, 0, height - 5);
        texture.needsUpdate = true;
    });
}


const loadButton = new qx.ui.form.Button("Load");

const updateStatus = function (message) {
  //  console.log(message)
    statusLabel.setValue(message);
};

const displayStatus = updateStatus;//_.throttle(updateStatus, 1000);

loadButton.addListener("execute", async function( e ) {

    loadButton.setEnabled(false);
    displayStatus("gathering volumes...");
    var suffix = ".nii.gz";
    await desk.FileSystem.traverseAsync( dataDirectory.getValue(),
        function (fileName) {
            if ((desk.FileSystem.getFileExtension(fileName) === "mhd") ||
                    (fileName.indexOf(suffix, fileName.length - suffix.length) !== -1)) {
                volumes.push (fileName);
                displayStatus("gathering volumes... " + volumes.length + " found ");
            }
        }
    );

    if ( cleanDB )  {

        var counter = 0;
        volumes = await Promise.all ( volumes.map( async function ( vol ) {
            displayStatus( "testing volume " + counter++ + "/" + volumes.length);
            if (vol.indexOf(suffix, vol.length - suffix.length) !== -1) {
                return vol;
            }
            var result = await desk.FileSystem.readFileAsync( vol );
            var test;
            try {
                var volume = MHD.parse( result );
                var index = 0;
                test = _.every( volume.DimSize, function ( size ) {
                    var length = index === 2 ? 0 : size * volume.ElementSpacing[ index ];
                    index++;
                    return size > 80 && length < 1000;
                });
            } catch (err) {
                return false;
            }
            return test ? vol : false;
    
        } ) );

        volumes = volumes.filter( v => v );

    }

    displayStatus(volumes.length + " valid  volumes");
    volumes = volumes.slice(0, dataSetSize);
//      volumes = _.sample(volumes, dataSetSize);

    volumes.sort();
    var massName = "data/volumes-" + desk.FileSystem.getFileName(dataDirectory.getValue()) +
        "-" + volumes.length + ".json";

    desk.FileSystem.writeFile(massName, JSON.stringify(volumes));
    loadVols();

});

async function loadVols() {

    await loadAll();
    pointsFile = [];
    matches = volumes.map( () => [] );
	displayStatus("Loading matches...");
	const cacheFile = "cache/" + getHash() + ".json";

	try {

        if ( !useCache ) throw ( 'no cache' );
        if ( !( await desk.FileSystem.existsAsync( cacheFile ) ) ) throw ( 'no cache' );
        const txt = await desk.FileSystem.readFileAsync( cacheFile );
        const response = JSON.parse( txt );
        console.log( "Use cached version of pairing and notes : " + cacheFile );
        matches = response.matches;
        pointsFile = response.pointsFile;

	} catch ( e ) {

        console.log("Generated pairing and notes");
        startTime = performance.now();
    	displayStatus("Computing matches...");
        for ( let volume of volumes ) await surf3d( volume );
        console.log("Save cache ...");
        await desk.FileSystem.writeJSONAsync ("cache/timings" + volumes.length + ".json", timings);
        await desk.FileSystem.writeJSONAsync ( cacheFile, { volumes, matches, pointsFile } );
        console.log( "FileSaved : " + cacheFile );

    }

    await solverWorker.setInput(matches, opts.edgeRemovalRatio, volumes);
    statusLabel.setValue("Ready to register all volumes. Press the \"cleanup\" button");
    //removeEdgesButton.setEnabled(true);
    resetRegistrationButton.setEnabled(true);
    registerButton.setEnabled(true);
    resetPositions();

}

const  viewer = new desk.THREE.Viewer( {
    orthographic : opts.orthographic,
    cameraFront : [ 0, 1, 0 ],
    cameraUp : [ 0, 0, 1 ]
//    cameraFront : [ 0, 0.94, -0.342 ],
//    cameraUp : [ 0, -0.342, 0.94 ]
});

if ( desk.auto ) viewer.getOptionsButton().setVisibility("excluded");
viewer.setPickMode(false);

const statusLabel = new qx.ui.basic.Label("");
const font = new qx.bom.Font( 30 );
statusLabel.setTextColor("blue");
statusLabel.setFont( font );
viewer.add( statusLabel, { left : "20%", top : 55 } );

viewer.addListener( "pick", onPick );
const resetView = _.throttle( viewer.resetView.bind( viewer ), 500 );
const cubeGeometry = new THREE.BoxGeometry( 1, 1, 1 );

function onPick(event) {

    const meshesSet = new Set();
    for ( let mesh of meshes ) meshesSet.add( mesh );
    viewer.setPickMode(false);
    pickButton.setEnabled(true);
    let currmesh = event.getData().object;
    while ( !meshesSet.has( currmesh ) ) currmesh = currmesh.parent;
    const meshid = meshes.indexOf( currmesh );
    currmesh.position.copy(getPosition(meshid));
    let pos = [];

    // get maximum number of matches
    let maxInliers = 0;
    meshes.forEach( function ( mesh, id ) {

        if ( meshid == id ) return;
        const transform = getMatrice(meshid, id);
        if (transform.fail === true) return;
        if (transform.inliers > maxInliers) maxInliers = transform.inliers;

    } );

    const bbox = new THREE.Box3();

    function setCubeToMesh( mesh, color ) {

        const cube = mesh.userData.cube;
        cube.material.color.set( color );
        cube.material.opacity = 0.5;
        mesh.remove( cube );
        bbox.setFromObject( mesh );
        mesh.add( cube );
        cube.visible = true;
        cube.scale.copy(bbox.max).sub(bbox.min);
        cube.position.copy(bbox.max).add(bbox.min).multiplyScalar(0.5).sub(mesh.position);

    }

    for ( let [ id, mesh ] of meshes.entries() ) {

        getPosition( id ).toArray( pos );
        let cube = mesh.userData.cube;

        if ( !cube ) {

            const cubeMaterial = new THREE.MeshLambertMaterial({ color: 0xFFFFFF });
            mesh.userData.cube = cube = new THREE.Mesh(cubeGeometry, cubeMaterial);
            cube.material.opacity = 0.5;
            cube.material.transparent = true;
            cube.visible = true;
            cube.renderOrder = 0;
            mesh.add( cube );

        }

        if ( meshid == id ) {
            setCubeToMesh( mesh, "yellow" );
            continue;
        }

        let transform = getMatrice( id, meshid );

        if ( ( id != meshid ) && transform && transform.fail) {
            mesh.position.fromArray( pos );
            cube.visible = true;
            setCubeToMesh( mesh, "black" );
            continue;
        }

        let scale = 1;
        if (meshid > id) {
            pos = numeric.add( pos, transform.translation );
            scale = transform.scale;
        }
        if (meshid < id) {
            pos = numeric.sub( pos, transform.translation );
            scale = 1 / transform.scale;
        }

        mesh.position.fromArray( pos );
        mesh.scale.setScalar( scale );

        if ( !transform || !transform.bboxA ) {
            console.log( "no transform");
            continue;
        }

        const bb = id < meshid ? transform.bboxA : transform.bboxB;
        const size = numeric.sub(bb.max, bb.min);
        const origin = bb.min;
        cube.visible = true;
        cube.scale.fromArray(size);
        cube.position.set(origin[0]+size[0]/2, origin[1]+size[1]/2, origin[2]+size[2]/2);

        let ratio = Math.pow( transform.inliers / maxInliers, 0.2 );
        cube.material.color.setRGB(ratio, 0, 1 - ratio);

    }

    viewer.render();

}

const removeEdgesButton = new qx.ui.form.ToggleButton( "cleanup" );
removeEdgesButton.addListener("execute", async function() {
    await startRegistration();
    removeEdgesButton.setEnabled(true);
});

const resetRegistrationButton = new qx.ui.form.Button("reset");
resetRegistrationButton.setEnabled(false);
resetRegistrationButton.addListener("execute", function() {

    statusLabel.setValue("Ready to register all volumes. Press the \"cleanup\" button");
    resetPositions();
    solverWorker.resetEdges();

} );

const registerButton = new qx.ui.form.Button("register");
registerButton.setEnabled(false);

registerButton.addListener("execute", function() {

    solverWorker.iterate(false, state => {
        updatePositions( state.positions );
        displayStatus( state.numberOfEdges + " edges" );
    } );
    
});


const pickButton = new qx.ui.form.Button("pick");
pickButton.addListener("execute", function() {
    viewer.setPickMode(true);
    pickButton.setEnabled(false);
    for ( let mesh of meshes ) {
        if ( !mesh.userData.cube ) continue;
        mesh.userData.cube.visible = false;
    }
    viewer.render();
} );


const logPosButton = new qx.ui.form.Button("log Position");
logPosButton.addListener("execute", function() {
    console.log(positions);
    console.log( positions.map( p => p.join( ", " ) ).join( "\n") );
} );

const unpickButton = new qx.ui.form.Button("reset pick");
unpickButton.addListener("execute", function() {
    for ( let mesh of meshes ) {
        mesh.visible = true;
        if ( mesh.userData.cube ) mesh.userData.cube.visible = false;
    }

    resetPositions();
});

// release memory
viewer.getWindow().addListener('close', function () {

    viewer.dispose();
    volumes.length = 0;
    matches.forEach(function (line) {
        line.length = 0;
    });
    matches.length = 0;

    solverWorker.terminate();
    viewer.getWindow().destroy();

});

viewer.addListener("drop", function(e) {
	if (e.supportsType("fileBrowser")) {
		loadFile(e.getData("fileBrowser").getSelectedFiles()[0]);
	}
});

async function loadFile( file ) {

	files = JSON.parse( await desk.FileSystem.readFileAsync( file ) );
    if ( !Array.isArray( files ) ) return;
    volumes = files.map( function ( file, index ) {
        if ( file.file ) file = file.file;
        return file;
    });

    if ( volumes.length < 10 ) {
        textureSamplingRatio = 1;
    } else {
        textureSamplingRatio = 0.2;
    }

    numberOfVolumes.setValue( files.length );
    loadVols();
}

const loadContainer = new qx.ui.container.Composite();
loadContainer.setLayout( new qx.ui.layout.HBox() );
const dataDirectory = new desk.FileField( dataDir );
dataDirectory.setWidth( 300 );
loadContainer.add( dataDirectory );

const numberOfVolumes = new qx.ui.form.Spinner();
numberOfVolumes.set({maximum: 3000, minimum: 1, value : dataSetSize});
numberOfVolumes.addListener("changeValue", function () {
    dataSetSize = numberOfVolumes.getValue();
    if ( dataSetSize > 10 ) width = Math.round(Math.sqrt(dataSetSize));
    else width = dataSetSize;
});

loadContainer.add( numberOfVolumes );
loadContainer.add( loadButton );
if ( !desk.auto && !loadDefaultFile ) viewer.add( loadContainer, { top : 10, right : 5 } );

var showSprites = new qx.ui.form.CheckBox('show valence');
showSprites.addListener("changeValue", displayAll);

const travellingButton = new qx.ui.form.Button("traveling");

travellingButton.addListener( "execute", async function () {

    const state1 = JSON.parse('{"zoomStart":0.10564663023679417,"zoomEnd":0.10564663023679417,"panStart":{"x":0.7067395264116576,"y":0.100789313904068},"panEnd":{"x":0.7067395264116576,"y":0.100789313904068},"eye":{"x":-212.35794042541238,"y":-6491.112456341291,"z":1508.5256828766226},"objectUp":{"x":0.011594266853438612,"y":0.2266986059974176,"z":0.9732159267731296},"objectPosition":{"x":-3510.8743948957235,"y":-1581.951260471882,"z":-394.6149375710809},"target":{"x":-3298.516454470311,"y":4909.161195869409,"z":-1903.1406204477034}}');
    const state2 = JSON.parse('{"zoomStart":0.10564663023679417,"zoomEnd":0.10564663023679417,"panStart":{"x":-0.071645415907711,"y":0.12082574377656345},"panEnd":{"x":-0.071645415907711,"y":0.12082574377656345},"eye":{"x":-212.35794042541056,"y":-6491.112456341291,"z":1508.5256828766226},"objectUp":{"x":0.011594266853438612,"y":0.2266986059974176,"z":0.9732159267731296},"objectPosition":{"x":13004.501386378697,"y":-2060.979345256498,"z":-132.26667752117032},"target":{"x":13216.859326804108,"y":4430.133111084793,"z":-1640.792360397793}}');
    const duration = 6;
    const fps = 30;
    const numberOfFrames = duration * fps;

    let frame = 0;
    viewer.getControls().setState(state1);
    viewer.render(true);
    if (snapshots) await viewer.snapshotAsync( {
        ratio : 3, path : "data/travelling/frame" + ( 10000 + frame ) + ".png" 
    } );

    while ( frame < numberOfFrames ) {
        frame++;
        viewer.getControls().interpolateState(state1, state2, frame / numberOfFrames);
        viewer.render(true);
        if (snapshots) {
            await viewer.snapshotAsync( {
                ratio : 3,
                path : "data/travelling/frame" + (10000 + frame) + ".png",
            });
        }

        await new Promise( res => setTimeout( res, 50 ) );
    }

} );

if ( !desk.auto ) viewer.add(travellingButton, {left : 0, bottom : 0});

if ( desk.auto ) {

    viewer.fillScreen();
	viewer.add(new qx.ui.embed.Html('Hubless 3D medical image bundle registration. More info : <a href="https://www.creatis.insa-lyon.fr/~valette/public/publication/agier-visapp-2016/">https://www.creatis.insa-lyon.fr/~valette/public/publication/agier-visapp-2016/</a>'), {top : 10, left : 10, width : 500});
	viewer.add(new qx.ui.embed.Html('Demo made with DESK. More info : <a href="https://www.creatis.insa-lyon.fr/~valette/public/project/desk">https://www.creatis.insa-lyon.fr/site/desk</a>'), {top : 30, left : 10, width : 500});

}

const buttonsContainer = new qx.ui.container.Composite();
buttonsContainer.setLayout( new qx.ui.layout.VBox() );
viewer.add( buttonsContainer, { bottom : 5, right : 5 } );
if ( !desk.auto ) buttonsContainer.add( logPosButton );
buttonsContainer.add( new qx.ui.basic.Label( "Registration") );
buttonsContainer.add( registerButton );
buttonsContainer.add( removeEdgesButton );
buttonsContainer.add( resetRegistrationButton );
if ( !desk.auto ) buttonsContainer.add( showSprites );
buttonsContainer.add( new qx.ui.basic.Label( "Picking") );
buttonsContainer.add( pickButton );
buttonsContainer.add( unpickButton );

viewer.add( new qx.ui.basic.Label( "contrast" ), { top : 80, right : 0 } );
const slider = new qx.ui.form.Slider( "vertical" );
slider.setHeight( 300 );
viewer.add( slider, { right : 10, top : 100 } );
slider.setValue( 100 - opts.contrast / 0.04 );
slider.addListener( "changeValue", () => {
    viewer.getScene().traverse( object => {
        const slice = object.userData?.viewerProperties?.volumeSlice;
        if ( !slice ) return;
        slice.setContrast( ( 100 - slider.getValue() ) * 0.04);
    } );
});

if ( desk.auto || loadDefaultFile ) loadFile( defaultFile );
