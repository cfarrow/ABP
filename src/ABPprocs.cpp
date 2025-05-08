#include <algorithm>
#include <iostream>
#include <queue>
#include <random>
#include <vector>

#include "Lattice.h" 
#include "Hexagonal.h"
#include "Square.h"
#include "Triangular.h"
#include "UJack.h"
#include "Cubic.h"
#include "BCC.h"
#include "FixedZ.h"
#include "SquRand.h"
#include "TriRand.h"
#include "CubRand.h"
#include "SWNetwork.h"
#include "Cubic4d.h"
#include "ABPprocs.h"
#include "rand.h"


using namespace std;

void collectToQueue (Lattice *L, queue< size_t > &Q, size_t cull) {

  for (size_t i = 0; i < L->getNumSites(); i++) {
      if (L->isActive(i) && L->getNumActiveNeighbors(i) < cull) {
		 Q.push(i);
		 L->setActiveLevel(i,0);
	 }
  }
  return;
}

/* ---------------------------------------------------------------*\

  Starting with the sites in the queue, make inactive all
  sites in the lattice that don't meet the culling condition.
  Returns the number of sites culled.

  Note that the Queue must already contain at least one site for 
  this alorithm to work. If a routine doesn't automatically put 
  sites in the culling queue, call the collectToQueue function 
  first.  

\* ---------------------------------------------------------------*/
size_t cullSites(Lattice *L, queue< size_t > &Q, size_t cull) {

  size_t culled = 0;
  size_t s1, s2;

  while (!Q.empty()) {		/* while queue still has sites  */

	  s1 = Q.front();
	  Q.pop();
	  culled++;			/* Keep track of how many sites have been culled.        */

      /* Now check all neighbor sites of the recently culled site and send them to the
         culling queue if they don't meet the culling condition.                         */
      for (size_t j = 0; j < L->getNumNeighbors(s1); j++) { 
		 s2 = L->getNbr(s1,j);
		 if ( L->isActive(s2) && L->getNumActiveNeighbors(s2) < cull) {
			 Q.push(s2);
			 L->setActiveLevel(s2,0);
		 }/* if */
	 }/* for j */
    
  }/* while */

  if (culled > L->getNumSites()) {
      cerr << "CullSites: Error! Culled more sites than available: " << culled << '.' << endl;
      exit (1);
    }

  return culled;
}

size_t removeActiveSite( Lattice *L, queue< size_t > &Q, size_t &t, size_t &k, bounded_rng_type& rng) {
	size_t j;

	/* burn sites (set absent) until an active site comes up */
	do {
		/* pick a random site */
		j = rng();

        /* If it is present, increment k and burn the site. If it is not
         * present, then it is not active. */
		if( L->isPresent(j) ) {
			k++;
			L->setPresentLevel(j, 0);
		}
	} while( !L->isActive(j) );

	/* An active site is found. It has already been set absent. Set inactive and push it onto the queue. */
	t++;
	L->setActiveLevel(j,0);
	Q.push(j);
	return j;
}

/* ---------------------------------------------------------------*\
   
   Various information-based functions.
   
\* ---------------------------------------------------------------*/

void promptLength( size_t &len ) {

	cout << "Enter the (even) side length of the lattice > ";
	cin >> len;

}

void promptCullParam( size_t &cull) {

	cout << "Enter the culling parameter > ";
	cin >> cull;

}

void promptNumSamples( size_t &num_samples) {

	cout << "How many samples would you like to run? (enter an integer) > ";
	cin >> num_samples; 
}


void GraphInfo::promptGraphType() {

	cout << "Enter the type of the lattice...\n";
	cout << HEX_LATTICE  << " = hexagonal" << endl;
	cout << SQU_LATTICE  << " = square" << endl;
	cout << TRI_LATTICE  << " = triangular" << endl;
	cout << JAK_LATTICE  << " = union jack" << endl;
	cout << CUB_LATTICE  << " = cubic" << endl;
	cout << BCC_LATTICE  << " = bcc" << endl;
	cout << C4D_LATTICE  << " = cubic4d" << endl;
	cout << FIZ_LATTICE  << " = random with fixed coordination" << endl;
	cout << SQR_LATTICE  << " = square/random" << endl;
	cout << TRR_LATTICE  << " = triangular/random" << endl;
	cout << CUR_LATTICE  << " = cubic/random" << endl;
	cout << SMW_LATTICE  << " = small world" << endl;
    cout << "> ";
	cin >> type;

    switch( type ){
        case HEX_LATTICE:
        case SQU_LATTICE:
        case TRI_LATTICE:
        case JAK_LATTICE:
            cout << "Enter the number of periodic directions (0-2) > ";
            cin >> pbc;
            if(pbc > 2) pbc = 2;
            if(pbc < 0) pbc = 0;
            break;
        case CUB_LATTICE:
        case BCC_LATTICE:
            cout << "Enter the number of periodic directions (0-3) > ";
            cin >> pbc;
            if(pbc > 3) pbc = 3;
            if(pbc < 0) pbc = 0;
            break;
        case C4D_LATTICE:
            cout << "Enter the number of periodic directions (0-4) > ";
            cin >> pbc;
            if(pbc > 4) pbc = 4;
            if(pbc < 0) pbc = 0;
            break;
        case FIZ_LATTICE:
            cout << "Enter the coordination of the lattice > ";
            cin >> coord;
            break;
        case SQR_LATTICE:
        case TRR_LATTICE:
        case CUR_LATTICE:
            cout << "Enter the fraction of random bonds > ";
            cin >> fb;
            break;
        case SMW_LATTICE:
            cout << "Enter the dimension of the network > ";
            cin >> dim;
            cout << "Enter the power-law exponent (alpha) > ";
            cin >> alpha;
            cout << "Enter the average coordination number (gamma) > ";
            cin >> gamma;
            break;
    }
}

/* output prompts */
void GraphInfo::printType( ostream &ostr ) {
	
	switch( type ) {
		case HEX_LATTICE:
			ostr << "hex";
			break;
		case SQU_LATTICE:
			ostr << "squ";
			break;
		case TRI_LATTICE:
			ostr << "tri";
			break;
		case JAK_LATTICE:
			ostr << "jak";
			break;
		case CUB_LATTICE:
			ostr << "cub";
			break;
		case BCC_LATTICE:
			ostr << "bcc";
			break;
		case C4D_LATTICE:
			ostr << "c4d";
			break;
		case FIZ_LATTICE:
			ostr << "fiz" << coord;
			break;
		case SQR_LATTICE:
			ostr << "sqr" << fb;
			break;
		case TRR_LATTICE:
			ostr << "trr" << fb;
			break;
		case CUR_LATTICE:
			ostr << "cur" << fb;
			break;
        case SMW_LATTICE:
            ostr << "smw_" << dim << "d_" << alpha << "a_" << gamma << "g";
            break;
        default:
            ostr << "err";
	}

    switch( type ) {
        case HEX_LATTICE:
        case SQU_LATTICE:
        case TRI_LATTICE:
        case JAK_LATTICE:
        case CUB_LATTICE:
        case BCC_LATTICE:
        case C4D_LATTICE:
            ostr << "_" << pbc << "pbc";
            break;
    }

}

void GraphInfo::printGraphInfo( ostream &ostr, size_t &nsites, size_t &mcull, size_t &num_samples) {

    // Print the graph type
    ostr << "# ";
    printType(ostr);
    ostr << "\n";
    // Length
    ostr << "# nsites = " << nsites << "\n";
    ostr << "# mcull = " << mcull << "\n";
    ostr << "# samples = " << num_samples << "\n";
}

Lattice* GraphInfo::getLattice(size_t len) {
	
	Lattice* Lat;
	switch( type ) {
		case HEX_LATTICE:
			switch(pbc){
				case 0: Lat = new Hexagonal<0>(len, type); break;
				case 1: Lat = new Hexagonal<1>(len, type); break;
				case 2: Lat = new Hexagonal<2>(len, type); break;
			}
			break;
		case SQU_LATTICE:
			switch(pbc){
				case 0: Lat = new Square<0>(len, type); break;
				case 1: Lat = new Square<1>(len, type); break;
				case 2: Lat = new Square<2>(len, type); break;
			}
			break;
		case TRI_LATTICE:
			switch(pbc){
				case 0: Lat = new Triangular<0>(len, type); break;
				case 1: Lat = new Triangular<1>(len, type); break;
				case 2: Lat = new Triangular<2>(len, type); break;
			}
			break;
		case JAK_LATTICE:
			switch(pbc){
				case 0: Lat = new UJack<0>(len, type); break;
				case 1: Lat = new UJack<1>(len, type); break;
				case 2: Lat = new UJack<2>(len, type); break;
			}
			break;
		case CUB_LATTICE:
			switch(pbc){
				case 0: Lat = new Cubic<0>(len, type); break;
				case 1: Lat = new Cubic<1>(len, type); break;
				case 2: Lat = new Cubic<2>(len, type); break;
				case 3: Lat = new Cubic<3>(len, type); break;
			}
			break;
		case BCC_LATTICE:
			switch(pbc){
				case 0: Lat = new BCC<0>(len, type); break;
				case 1: Lat = new BCC<1>(len, type); break;
				case 2: Lat = new BCC<2>(len, type); break;
				case 3: Lat = new BCC<3>(len, type); break;
			}
			break;
		case C4D_LATTICE:
			switch(pbc){
				case 0: Lat = new Cubic4d<0>(len, type); break;
				case 1: Lat = new Cubic4d<1>(len, type); break;
				case 2: Lat = new Cubic4d<2>(len, type); break;
				case 3: Lat = new Cubic4d<3>(len, type); break;
				case 4: Lat = new Cubic4d<4>(len, type); break;
			}
			break;
		case FIZ_LATTICE:
			return new FixedZ(len, type, coord);
			break;
		case SQR_LATTICE:
			return new SquRand(len, type, fb);
			break;
		case TRR_LATTICE:
			return new TriRand(len, type, fb);
			break;
		case CUR_LATTICE:
			return new CubRand(len, type, fb);
			break;
		case SMW_LATTICE:
			return new SWNetwork(len, type, dim, alpha, gamma);
			break;
        default:
            return NULL;
	}
	return Lat;
}

/* Process command line arguments 
 *
 * If the command line cannot be processed, the user will be prompted for info.
 *
 * the order in which information is passed:
 * type
 * len/nsites
 * mcull
 * num_samples
 * pbc
 * coord
 * fb
 * dim
 * alpha
 * gamma
 *
 * FIXME - make this work with command line switches
 *
 */
void processCommandLine(int argc, char **argv, size_t &len, size_t &mcull, size_t &num_samples, GraphInfo &gi)
{

	if( argc >= 5 ) { 
		gi.type = atoi(argv[1]);
		len = atoi(argv[2]);
		mcull = atoi(argv[3]);
		num_samples = atoi(argv[4]);
	}
    if( argc >= 6 ) { gi.pbc = atoi(argv[5]); }
    if( argc >= 7 ) { gi.coord = atoi(argv[6]); }
    if( argc >= 8 ) { gi.fb = atof(argv[7]); }
    if( argc >= 9 ) { gi.dim = atoi(argv[8]); }
    if( argc >= 10) { gi.alpha = atof(argv[9]); }
    if( argc == 11 ) { gi.gamma = atof(argv[10]); }
    if( argc > 11 || argc < 5 ) {
		gi.promptGraphType();
        promptLength( len );
        promptCullParam( mcull );
        promptNumSamples( num_samples );
        cout << "\nThanks!" << endl;
    }
}