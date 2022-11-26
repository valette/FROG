
var bbox = new THREE.Box3();

THREE.SpriteHelper = function ( canvas ) {
    texture = new THREE.Texture( canvas );
    texture.minFilter = THREE.LinearFilter; // to avoid a warning on texture size
    texture.needsUpdate = true;
    var material = new THREE.SpriteMaterial( { map :texture, depthTest : false } );
    var sprite = new THREE.Sprite( material );
    return sprite;
};

THREE.TextSpriteHelper = function (obj, text, options) {
    options = options || {};
    var texture, canvas, context;
    var sprite = obj.userData.sprite;
    var size = options.size || 100;
    var height = 64;
    text = text.trim();
    if (!sprite) {
        canvas = document.createElement('canvas');
        canvas.height = height;
        context = canvas.getContext("2d");
        context.font = height + 'px bold Arial';
        var width = context.measureText( text + '' ).width;
        canvas.width =  width;
        sprite = new THREE.SpriteHelper( canvas );
        bbox.setFromObject( obj );
        if ( bbox.isEmpty() ) bbox.expandByPoint( obj.position );
        bbox.translate( obj.position.clone().negate() );
        sprite.position.copy( bbox.getCenter( new THREE.Vector3() ) );
        sprite.position.z = bbox.max.z;
        sprite.scale.x = size * width / height;
        sprite.scale.y = size;
        obj.add(sprite);
    }

    context = canvas.getContext( '2d' );
    texture = sprite.material.map;
    canvas = texture.image;
    context.clearRect(0, 0, canvas.width, canvas.height);
    context.fillStyle = options.backgroundColor || 'rgba(0, 0, 0, 0.5)';
    context.fillRect(0, 0, canvas.width, canvas.height);
    context.fillStyle = options.color || 'yellow';
    context.font = height + 'px bold Arial';
    context.fillText( text, 0, height - 7 );
    if (options.strokeColor) {
        context.strokeStyle = options.strokeColor;
        context.lineWidth = 2;
        context.strokeText( text, 0, height - 7 );
    }
    texture.needsUpdate = true;
    return sprite;
};

// getColorFromValue is a function which should return an rgb color
// where each color component is between 0 and 1
THREE.ColorBarSpriteHelper = function ( getColorFromValue, options ) {
    options = options || {};
    var scale = options.scale || 10;
    var fontHeight = options.fontHeight || 100;
    var width = options.width || 50;
    var height = options.height || 500;
    var yOffset = options.yOffset || 10;
    var max = options.max || 1;

    var bCanvas = document.createElement('canvas');
    var bHeight = height + fontHeight;
    bCanvas.width = 400;
    bCanvas.height = bHeight;
    var bContext = bCanvas.getContext('2d');
    bContext.font = 100 + 'px bold Arial';
    bContext.fillText( '0', width + 10, bHeight  - yOffset);
    bContext.fillText( max.toPrecision(3) , width + 10, fontHeight  - yOffset);

    var canvas = document.createElement('canvas');
    canvas.width = width;
    canvas.height = height;
    var context = canvas.getContext("2d");
    var imageData = context.createImageData(width, height);
    var data = imageData.data;
    for (var i = 0; i < height; i++) {
        for (var j = 0; j < width; j ++) {
            var col = getColorFromValue(1 - i / height);
            var offset = 4 * i * width + 4 * j;
            data[ offset ] = 255 * col[0];
            data[ offset + 1 ] = 255 * col[1];
            data[ offset + 2 ] = 255 * col[2];
            data[ offset + 3 ] = 255;
        }
    }
    bContext.putImageData(imageData, 0, fontHeight / 2);
    var sprite = new THREE.SpriteHelper( bCanvas );
    sprite.scale.multiplyScalar(bHeight);
    sprite.frustumCulled = false;
    return sprite;
};
