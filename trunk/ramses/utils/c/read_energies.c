// read_energies
//=============================================================================
// Author: Michael Rieder (2015) rieder@physik.uzh.ch
//
// parse log file for energy outputs
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
	// Main step
	const char * mainstep_string = "Main step=";
	char * mainstep_pos;
	// econs
	const char * econs_string = "econs=";
	char * econs_pos;
	double	econs;
	// epot
	const char * epot_string = "epot=";
	char * epot_pos;
	double	epot;
	// ekin
	const char * ekin_string = "ekin=";
	char * ekin_pos;
	double	ekin;
	// eint
	const char * eint_string = "eint=";
	char * eint_pos;
	double	eint;
	// emag
	const char * emag_string = "emag=";
	char * emag_pos;
	double	emag;
	// t
	const char * time_string = "t=";
	char * time_pos;
	double	time;
	// dt
	const char * deltat_string = "dt=";
	char * deltat_pos;
	double	deltat;
	// a
	const char * aexp_string = "a=";
	char * aexp_pos;
	double	aexp;

	// read new line
	while( fgets( line_buffer, max_line_length, logfile ) ) {
		mainstep_pos = strstr( line_buffer, mainstep_string );
		if ( mainstep_pos ) {
			// econs
			econs_pos = strstr( line_buffer, econs_string );
			econs = atof( econs_pos + strlen(econs_string) );
			// epot
			epot_pos = strstr( line_buffer, epot_string );
			epot = atof( epot_pos + strlen(epot_string) );
			// ekin
			ekin_pos = strstr( line_buffer, ekin_string );
			ekin = atof( ekin_pos + strlen(ekin_string) );
			// eint
			eint_pos = strstr( line_buffer, eint_string );
			eint = atof( eint_pos + strlen(eint_string) );

			// read next line
			fgets( line_buffer, max_line_length, logfile );
			// try emag
			emag_pos = strstr( line_buffer, emag_string );
			if ( emag_pos ) {
				emag = atof( emag_pos + strlen(emag_string) );
				// read next line (time and aexp)
				fgets( line_buffer, max_line_length, logfile );
			}
			else emag = 0.0;

			// time
			time_pos = strstr( line_buffer, time_string );
			time = atof( time_pos + strlen(time_string) );
			// deltat
			deltat_pos = strstr( line_buffer, deltat_string );
			deltat = atof( deltat_pos + strlen(deltat_string) );
			// aexp
			aexp_pos = strstr( line_buffer, aexp_string );
			aexp = atof( aexp_pos + strlen(aexp_string) );

			// output t - a - e values
			printf( "%f %f %e %e %e %e %e %e\n", time, aexp, deltat,
						econs, epot, ekin, eint, emag );
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

