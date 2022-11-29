"use strict";

var Heap = require ('heap');


var laplaceSolver = {
	worker : function () {
	    const heapScript = 'FROG/lib/heap.js';
	    desk.FileSystem.getFileURL( heapScript ); // statification hack
		return operative( laplaceSolver, [ heapScript] );
	},

	laplacianValues : null,
	values : null,
	newValues : null,
	laplacianWeights : null,

	solve : function() {
		var matches = this.matches;
		var i, j, i3, j3, weight, pos;

		if(!this.isGraphConnex(true)) {
			console.error("Graph non connexe");
		}

		if (!this.laplacianValues) {
			// init arrays
			this.laplacianValues = new Float64Array(matches.length * 4);
			this.values = new Float64Array(matches.length * 4);
			this.newValues = new Float64Array(matches.length * 4);
			this.laplacianWeights = new Uint16Array(matches.length);
		}

		for (i = 0 ; i < this.values.length ; i++) {
			this.values[i] = 0;
		}

		var laplacianValues = this.laplacianValues;
		var values = this.values;
		var newValues = this.newValues;
		var laplacianWeights = this.laplacianWeights;

		// reset arrays
		for (i = 0 ; i < laplacianValues.length ; i++) {
			laplacianValues[i] = 0;
		}

		for (i = 0 ; i < laplacianWeights.length ; i++) {
			laplacianWeights[i] = 0;
		}

		// compute laplacian
		for (i = 0 ; i < matches.length ; i++) {
			i3  = 4 * i;
			for (j = i + 1 ; j < matches.length ; j++) {
				if (matches[i][j].fail === true) {
					continue;
				}
				j3 = 4 * j;
				var translation = matches[i][j].translation;
				laplacianValues[i3] += translation[0];
				laplacianValues[i3 + 1] += translation[1];
				laplacianValues[i3 + 2] += translation[2];
				laplacianValues[i3 + 3] += Math.log(matches[i][j].scale);
				laplacianValues[j3] -= translation[0];
				laplacianValues[j3 + 1] -= translation[1];
				laplacianValues[j3 + 2] -= translation[2];
				laplacianValues[j3 + 3] -= Math.log(matches[i][j].scale);
				laplacianWeights[i] ++;
				laplacianWeights[j] ++;
			}
		}

		// normalize laplacian
		for (i = 0 ; i < matches.length ; i++) {
			weight = laplacianWeights[i];
			laplacianValues[4 * i] /= weight;
			laplacianValues[4 * i + 1] /= weight;
			laplacianValues[4 * i + 2] /= weight;
			laplacianValues[4 * i + 3] /= weight;
		}

		for (var loop = 0; loop < 10 * matches.length ; loop++) {
			// reset positions
			for (i = 0 ; i < 4 * matches.length ; i++) {
				newValues[i] = 0;
			}

			for (i = 0 ; i < matches.length ; i++) {
				i3  = 4 * i;
				for (j = i + 1 ; j < matches.length ; j++) {
					if (matches[i][j].fail === true) {
						continue;
					}

					j3 = 4 * j;
					newValues[i3] += values[j3];
					newValues[i3 + 1] += values[j3 + 1];
					newValues[i3 + 2] += values[j3 + 2];
					newValues[i3 + 3] += values[j3 + 3];

					newValues[j3] += values[i3];
					newValues[j3 + 1] += values[i3 + 1];
					newValues[j3 + 2] += values[i3 + 2];
					newValues[j3 + 3] += values[i3 + 3];
				}
			}

			for (i = 0 ; i < matches.length ; i++) {
				i3 = 4 * i;
				weight = laplacianWeights[i];
				newValues[i3] = laplacianValues[i3] + newValues[i3] / weight;
				newValues[i3 + 1] = laplacianValues[i3 + 1] + newValues[i3 + 1] / weight;
				newValues[i3 + 2] = laplacianValues[i3 + 2] + newValues[i3 + 2] / weight;
				newValues[i3 + 3] = laplacianValues[i3 + 3] + newValues[i3 + 3] / weight;
			}
			// swap arrays
			var temp = newValues;
			newValues = values;
			values = temp;
		}
		this.values = values;
		this.newValues = newValues;
		return;
	},

	// members for connexity trests
	tags : null,
	stack : null,

	isGraphConnex : function ( debug, fix ) {

		const matches = this.matches;

		if ( !this.tags ) {

			this.tags = new Uint8Array( matches.length );
			this.stack = new Uint16Array( matches.length * matches.length );

		}

		const tags = this.tags;
		const stack = this.stack;
		for ( let i = 0 ; i < matches.length ; i++) tags[ i ] = 0;
		let stackPointer = 1;
		stack[ 0 ] = 0;

		while (stackPointer--) {

			let s = stack[ stackPointer ];
			tags[ s ] = 1;

			for ( let t = 0 ; t < matches.length ; t++) {

				if ((t !=s ) && (matches[Math.min(t, s)][Math.max(t, s)] === undefined))
					console.log("missing couple ", t, s);

				if ( (t !== s && tags[t] === 0) &&
				    ( matches[Math.min(t, s)][Math.max(t, s)].fail !== true) )
    					stack[stackPointer++] = t;

			}

		}

		let connex = true;

		for ( let i = 0; i < tags.length; i++) {

			if ( tags[ i ] === 0 ) {

				if (debug) {

					console.log("volume " + i + " is disconnected ");
					console.log(this.volumes[i]);

				}

				connex = false;

			}

		}

        if ( !connex && fix ) {

            console.log( "fix disconnected volumes:" );
            const disconnected = [];
            let firstConnected = -1;

    		for ( let i = 0; i < tags.length; i++) {

    			if ( tags[ i ] === 0 ) disconnected.push( i );
    			else  if ( firstConnected < 0 ) firstConnected = i;

    		}

            console.log( disconnected );
            console.log( "First connected volume : " + firstConnected );

            for ( let i of disconnected ) {

                matches[ Math.min( firstConnected, i)][Math.max(firstConnected, i) ] = 
                {
                    fail : false,
                    scale : 1,
                    translation : [ 0, 0, 0 ]
                };

            }

        }

		return connex;

	},

	// members for edge removal
	heap : null,

	removeEdgeBatch : function () {

		if ( this.numberOfEdges === 0 ) this.getNumberOfEdges();
		const matches = this.matches;

		if ( !this.heap ) {

			this.heap = new Heap(function(a, b) {
				if ( a.note !== b.note ) return a.note - b.note;
				if ( a.i !== b.i ) return a.i - b.i;
				return a.j - b.j;
			});

			for ( let i = 0 ; i < matches.length ; i++) {
				for ( let j = i + 1 ; j < matches.length ; j++ ) {
					if (matches[i][j].fail === true) {
						continue;
					}

					this.heap.push( { note : matches[ i ][ j ].inliers, i, j } );
				}
			}
		}

		const heap = this.heap;
		var toRemove = [];

		for ( let i = 0; i < this.edgeRemovalRatio * heap.size(); i++) {

			if ( heap.empty() ) break;
			toRemove.push(heap.pop());

		}

		function removeEdge( edge ) {
			matches[ edge.i ][ edge.j ].fail = true;
		}

		function addEdge( edge ) {
			const e = matches[edge.i][edge.j];
			e.fail = e.f;
		}

		const stack = [ toRemove ];

		while ( stack.length ) {

			const edges = stack.pop();
			for( let edge of edges ) removeEdge( edge );

			if ( !this.isGraphConnex() ) {

    			for( let edge of edges ) addEdge( edge );
				if ( edges.length === 1 ) continue;
				const half = Math.round( edges.length / 2 );
				stack.push( edges.slice( half, edges.length ) );
				stack.push( edges.slice( 0, half ) );

			}

		}

	},

	numberOfEdges : 0,
	valence : null,

	getNumberOfEdges : function () {

		const matches = this.matches;
		if ( !this.valence ) this.valence = new Uint32Array(matches.length);
		this.numberOfEdges = 0;
		for( let i = 0 ; i < matches.length; i++) this.valence[ i ] = 0;

		for( let i = 0 ; i < matches.length ; i++) {
			for (let j = i + 1 ; j < matches.length ; j++) {
				if ( !matches[i][j].fail ) {
					this.numberOfEdges++;
					this.valence[ i ]++;
					this.valence[ j ]++;
				}
			}
		}
	},

	iterate : function ( removeEdges, callback ) {
/*		if (!this.isGraphConnex(true)) {
			callback({error : "graph non connex"});
			return;
		}*/

        this.isGraphConnex(true, true );
		if (removeEdges) this.removeEdgeBatch();
		this.solve();
        // set first coordinates to 0;
        var anchor = [this.values[0], this.values[1], this.values[2], this.values[3]];

		// normalize laplacian
		for (let i = 0 ; i < this.matches.length ; i++) {
		    for (let j = 0; j < 4; j++) {
    			this.values[4 * i + j] -= anchor[j];
		    }
		}

		this.getNumberOfEdges();

		callback({positions : this.values, valences : this.valence, numberOfEdges : this.numberOfEdges});
	},

	setInput : function ( matches, edgeRemovalRatio, volumes, callback ) {
		this.matches = matches;
		this.laplacianValues = 0;
		this.edgeRemovalRatio = edgeRemovalRatio;
		this.volumes = volumes;
		for ( let line of this.matches ) {
			for ( let t of line ) {
				if (!t) continue;
				t.f = t.fail;
			}
		}
		if (typeof callback === "function") callback();
	},

	resetEdges : function ( callback ) {

		this.numberOfEdges = this.heap = 0;

		for ( let line of this.matches ) {
			for ( let t of line ) {
				if (!t) continue;
				t.fail = t.f;
			}
		}

		callback();

	}

};

if ( typeof define === 'function' && define.amd ) {

		define( 'laplaceSolver', laplaceSolver );

} else if ( 'undefined' !== typeof exports && 'undefined' !== typeof module ) {

		module.exports = laplaceSolver;

} else {
  self.laplaceSolver = laplaceSolver;
}
