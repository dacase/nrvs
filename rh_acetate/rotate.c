#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "nabcode.h"
extern char NAB_rsbuf[];
static int mytaskid, numtasks;

static MOLECULE_T *m;


static MATRIX_T mat;

static REAL_T phi;


int main( argc, argv )
	int	argc;
	char	*argv[];
{
	nabout = stdout; /*default*/

	mytaskid=0; numtasks=1;
m = getpdb( "dac3.pdb", NULL );





phi = torsion( m, "::O3", "::C3", "::C4", "::H4" );
NAB_matcpy( mat, rot4( m, "::C3", "::C4",  - phi ) );
transformmol( mat, m, "::H4,H5,H6" );

phi = torsion( m, "::O4", "::C3a", "::C4a", "::H4a" );
NAB_matcpy( mat, rot4( m, "::C3a", "::C4a",  - phi ) );
transformmol( mat, m, "::H4a,H5a,H6a" );

phi = torsion( m, "::O1", "::C1", "::C2", "::H1" );
NAB_matcpy( mat, rot4( m, "::C1", "::C2",  - phi ) );
transformmol( mat, m, "::H1,H2,H3" );

phi = torsion( m, "::O2", "::C1a", "::C2a", "::H1a" );
NAB_matcpy( mat, rot4( m, "::C1a", "::C2a",  - phi ) );
transformmol( mat, m, "::H1a,H2a,H3a" );

putpdb( "dac4.pdb", m, NULL );



	exit( 0 );
}
