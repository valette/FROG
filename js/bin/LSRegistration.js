#!/usr/bin/env node

'use strict';

const async          = require ( 'async' ),
	  fs             = require ( 'fs' ),
	  LSRegistration = require ( __dirname + '/../lib/LSRegistration.js' ),
	  { promisify }  = require ( 'util' );

( async function () {

try {

	const inputFile = process.argv[ 2 ];
	console.log( "Loading : " + inputFile );
	console.time( "total computation" );
	const obj = require( inputFile );
	const volumes = obj.volumes || obj;
	console.log( "Computing registration..." );
	const registration = new LSRegistration();
	registration.registerAsync = promisify( registration.register );
	const res = await registration.registerAsync( volumes );
	res.volumes = volumes;
	fs.writeFileSync( "registration.json" , JSON.stringify( res ) );
	console.log( "... Done" );
	console.timeEnd( "total computation" );
	process.exit( 0 );

} catch ( e ) {

	console.log( e );
	process.exit( 1 );

}

} )();
