// read_emag
//=============================================================================
// Author: Michael Rieder (2013) rieder@physik.uzh.ch
//
// parse log for emag
//=============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>


// global variables

const char *	input_filename = 0;
const int	max_line_length = 128;

// function declarations

void	read_logfile( FILE *logfile );
FILE*	open_logfile( const char *filename );
void	close_logfile( FILE *logfile );

void	parse_args( int argc, char *argv[] );


// function definitions

int main( int argc, char *argv[] )
{
	FILE *logfile;

	// check the arguments
	parse_args( argc, argv );	

	if ( !input_filename )
		return EXIT_FAILURE;

	// file access
	logfile = open_logfile( input_filename );

	if ( logfile ) {
		read_logfile( logfile );
		close_logfile( logfile );
	}

	return EXIT_SUCCESS;
}

void read_logfile( FILE *logfile )
{
	char line_buffer[max_line_length];
	// emag
	const char * emag_string = "emag=";
	char * emag_pos;
	double	emag;
	// t
	const char * time_string = "t=";
	char * time_pos;
	double	time;
	// a
	const char * aexp_string = "a=";
	char * aexp_pos;
	double	aexp;

	// read new line
	while( fgets( line_buffer, max_line_length, logfile ) ) {
		emag_pos = strstr( line_buffer, emag_string );
		if ( emag_pos ) {
			emag = atof( emag_pos + strlen(emag_string) );

			// read next line
			fgets( line_buffer, max_line_length, logfile );
			// time
			time_pos = strstr( line_buffer, time_string );
			time = atof( time_pos + strlen(time_string) );
			// aexp
			aexp_pos = strstr( line_buffer, aexp_string );
			aexp = atof( aexp_pos + strlen(aexp_string) );

			// output t - a - e values
			printf( "%f %f %e\n", time, aexp, emag );
		}
	}
}

FILE* open_logfile( const char *filename )
{
	FILE *logfile;

	logfile = fopen( input_filename, "r" );
	if ( !logfile )
		printf( "error opening %s\n", input_filename );

	return logfile;
}

void close_logfile( FILE *logfile )
{
	if ( logfile )
		fclose( logfile );
}

void parse_args( int argc, char *argv[] )
{
	int i;

	if ( argc < 2 ) {
		printf( "Usage: %s filename\n", argv[0] );
	}

	for (i=1;i<argc;i++) {
		if ( argv[i][0] == '-' ) {
			printf( "unrecognized command: %s\n", argv[i] );
		}
		else input_filename = argv[i];
	}
}

