/* -------------------------------------------------------------------------------*\ 
 09/23/04 Calculates avalanche size distribution for all values of t. Outputs
 		distribution in array format.
\* -------------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <cmath>
#include <queue>
#include <stack>
#include <vector>
#include <numeric>
#include "Lattice.h"  /* compile with "Lattice.cpp" */
#include "ABPprocs.h" /* compile with "ABPprocs.cpp" */
#include "MersenneTwister.h" /* compile with "ABPprocs.cpp" */

using namespace std;

int main(int argc, char **argv) {

	queue< size_t > CullingQ; /* The culling queue. */

	size_t type; /* The type of lattice */
	size_t len; /* The length of the lattice */
    size_t coord = 6;   /* coordination of fixed-coordination lattice */
    short pbc = 0; /* periodic boundary conditions */
    double fb = 0; /* Fraction of random bonds in TriRand lattice */

	size_t mcull; /* The culling parameter */
	size_t num_samples;
	size_t num_sites;  /* The number of sites in the lattice */

	if( argc >= 5 ) { 
		type = atoi(argv[1]);
		len = atoi(argv[2]);
		mcull = atoi(argv[3]);
		num_samples = atoi(argv[4]);
	}
    if( argc >= 6 ) { pbc = atoi(argv[5]); }
    if( argc >= 7 ) { coord = atoi(argv[6]); }
    if( argc == 8 ) { fb = atof(argv[7]); }
	if( argc > 8 || argc < 5 ) {
		cout << "This program performs the ABP process and calculates\n";
	    cout <<	"the cumulative avalanche size distribution for the\n";
		cout << "entire range of t." << endl;
		promptGraphType( type, coord, pbc, fb );
		promptLength( len );
		promptCullParam( mcull );
		promptNumSamples( num_samples );
		cout << "\nThanks!" << endl;
	}
	
    Lattice *Lat1 = getLattice( len, type, coord, pbc, fb);
	num_sites = Lat1->getNumSites();
	len = Lat1->getLength();

	/* seed and name the rng */
	MTRand();
	MTRand mtRNG;

	size_t a_size = 0;
	size_t a_max = 0;
	double a_max_avg = 0.0;
	double t_e = 0.0;
	double k_e = 0.0;
	double tau = 0.0;
	size_t t = 0;
	size_t t_max = 0;
	size_t k = 0;
	size_t a_range = (size_t) (num_sites);
	size_t t_range = 100;
	/* 1/t_range is the interval of the chosen t-values */
	
	/* dynamically allocate the distribution double-array */
	double *d_ref = new double [ a_range * t_range ];
	for( size_t i = 0; i < a_range * t_range; i++) {
		d_ref[i] = 0.0;
	}

	/* a_dist[t][a] gives the number of avalanches of size a at time t */
	double **a_dist = new double* [t_range];
	for( size_t i = 0; i < t_range; i++ ) {
		a_dist[i] = & d_ref[ a_range * i ];
	}	

	for( size_t smpl = 0; smpl < num_samples; smpl++ ) {
	
		cout << "sample " << smpl+1 << endl;
		Lat1->activateSites();

		t = 0;
		tau = 0.0;
		k = 0;

		/* Evacuate the lattice, recording the history avalanche sizes
		   in the vector avalanches. */

		while( Lat1->getNumActive() > 0 ) {
			removeActiveSite(Lat1, CullingQ, t, k, mtRNG);
			a_size = cullSites(Lat1, CullingQ, mcull) - 1;
			tau = (1.0*t)/num_sites;
			for( size_t j = 1; j < t_range; j++) {
				if( tau <= (1.0*j/t_range) ) {
					a_dist[j][a_size]++;
				}
			}
			if( a_size > a_max ) { a_max = a_size; }
		} /* while */
		a_max_avg += a_max;

		if( t > t_max ) { t_max = t; }
		t_e += t;
		k_e += k;

	} /* smpl */
	a_max_avg /= num_samples;

	/* rescale find t_range_max so no partial statistics are used */
	size_t t_range_max = static_cast<size_t>(t_range*t_max*1.0/num_sites);

	/* normalize the data to the number of samples and size of the lattice */
	for( size_t i = 1; i <= t_range_max; i++ ) {
		for( size_t j = 0; j < a_max_avg; j++ ) {
			a_dist[i][j] /= num_samples*num_sites;
			/* must divide a_dist[i][j] by its actual t-value */
			a_dist[i][j] /= num_sites * (1.0*i) / t_range;
		}
	}

	t_e /= num_samples;
	k_e /= num_samples;

	cout << "\nPrinting data to files.";

    stringstream filename;
    filename << "Ca_dist_t_" << mcull;
	printType( filename, type, coord, pbc, fb );
    filename << ".dat";
	ofstream adt_file( filename.str().c_str(), ios::app );
	adt_file << "#Cumulative Avalanche Numbers. ABPCDist_t\n" ;
	adt_file << "#Data presented in matrix format n[a][t] (columns are time).\n";
	adt_file << "# ";
	printType( filename, type, coord, pbc, fb );
	
	adt_file << "\n#m = " << mcull << " samples = " << num_samples << endl;
	adt_file << "# L = " << len << " t_e = " << t_e/num_sites << " k_e = " << k_e/num_sites << endl;

	adt_file << "#";
	for( size_t i = 0; i <= t_range_max; i++ ) {
		adt_file << (1.0*i) / t_range << '\t';
	}
	adt_file << endl;
	for( size_t j = 0; j < a_max_avg; j++ ) {

		for( size_t i = 0; i <= t_range_max; i++ ) {
			if( i == 0 ) { adt_file << j << '\t'; }
			else {
				adt_file << a_dist[i][j] << '\t' ;
			}
		}
		adt_file << endl;

	} 
	adt_file << endl;

	cout << endl;

	return 0;
}
