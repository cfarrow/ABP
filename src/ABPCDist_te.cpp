/* -------------------------------------------------------------------------------*\ 
 09/23/04   Calculates avalanche size distribution for all values of t. Outputs
            distribution in array format.
 03/05/07   Modified code to only record information for t_e.

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
#include "MersenneTwister.h"

using namespace std;

int main(int argc, char **argv) {

	queue< size_t > CullingQ; /* The culling queue. */

	size_t len; /* The length of the lattice */
	size_t mcull; /* The culling parameter */
	size_t num_samples;
	size_t num_sites;  /* The number of sites in the lattice */
    GraphInfo gi; /* Class for holding graph information */

    cout << "This program performs the ABP process and calculates the avalanche size \n";
    cout << "distribution over the entire evacuation process." << endl;
    processCommandLine(argc, argv, len, mcull, num_samples, gi);
    Lattice *Lat1 = gi.getLattice(len);
	num_sites = Lat1->getNumSites();
    len = Lat1->getLength();


	/* seed and name the rng */
	MTRand();
	MTRand mtRNG;

	size_t a_size = 0;
	size_t a_max = 0;
	double t_e = 0.0;
	double k_e = 0.0;
	size_t t = 0;
	size_t k = 0;
	size_t a_range = (size_t) (num_sites);
	
	double *a_dist = new double [a_range];
	for( size_t i = 0; i < a_range; i++ ) {
		a_dist[i] = 0.0;
	}	

	for( size_t smpl = 0; smpl < num_samples; smpl++ ) {
	
		cout << "sample " << smpl+1 << endl;
		Lat1->activateSites();

		t = 0;
		k = 0;

		/* Evacuate the lattice, recording the history avalanche sizes
		   in the vector a_dist. */

		while( Lat1->getNumActive() > 0 ) {
			removeActiveSite(Lat1, CullingQ, t, k, mtRNG);
			a_size = cullSites(Lat1, CullingQ, mcull) - 1;
            ++a_dist[a_size];
            if(a_size > a_max) a_max = a_size;
		} /* while */

		t_e += t;
		k_e += k;

	} /* smpl */
	t_e /= num_samples*num_sites;
	k_e /= num_samples*num_sites;

    /* Normalize the data so that distribution sums to 1 over avalanche size.
     */
    double atot = 0.0;
	for( size_t i = 0; i < a_range; i++ ) {
        atot += a_dist[i];
	}
	for( size_t i = 0; i < a_range; i++ ) {
        a_dist[i] /= atot;
	}

	cout << "\nPrinting data to files.";

    stringstream filename;
    filename << "Ca_dist_te_" << mcull;
	gi.printType(filename); /* Append graph info to file name */
    filename << ".dat";
	ofstream adt_file( filename.str().c_str(), ios::app );
	adt_file << "#Cumulative Avalanche Numbers. ABPCDist_te\n" ;
	adt_file << "# ";
    gi.printGraphInfo(adt_file, num_sites, mcull, num_samples);
	adt_file << "# t_e = " << t_e << " k_e = " << k_e << endl;

	adt_file << "#" << endl;
	for( size_t j = 0; j <= a_max; j++ ) {

        adt_file << j << '\t' << a_dist[j] << endl;

	} 
	adt_file << endl;

	cout << endl;

	return 0;
}
