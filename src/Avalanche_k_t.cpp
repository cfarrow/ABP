/* -------------------------------------------------------------------------------*\ 
 08/25/04 Finds avalanche size vs k and t also calculates k_vs_t.
 11/01/04 Includes "negative avalanches in k".
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
#include "rand.h"

using namespace std;

int main(int argc, char **argv) {

	queue< size_t > CullingQ; /* The culling queue. */

	size_t len; /* The length of the lattice */
	size_t mcull; /* The culling parameter */
	size_t num_samples;
	size_t num_sites;  /* The number of sites in the lattice */
    GraphInfo gi; /* Class for holding graph information */

    cout << "This program performs the ABP process and calculates the avalanche size as a\n";
    cout << "function of t and k for a given lattice type." << endl;
    processCommandLine(argc, argv, len, mcull, num_samples, gi);
    Lattice *Lat1 = gi.getLattice(len);
	num_sites = Lat1->getNumSites();
    len = Lat1->getLength();

	/* seed and name the rng */
	bounded_rng_type rng = makeRNG(7, 0, num_sites-1);

	size_t a_size = 0;
	double t_e = 0.0;
	double k_e = 0.0;
	size_t k_max = 0;
	size_t k_last = 0;
	size_t t_max = 0;
	size_t t = 0;
	size_t k = 0;
	vector< double > avalanches_t;
	vector< double > avalanches_k;
	vector< double > k_vs_t;

	avalanches_t.resize( num_sites+1 );
	avalanches_k.resize( num_sites+1 );
	k_vs_t.resize( num_sites+1 );
	fill( avalanches_t.begin(), avalanches_t.end(), 0.0 );
	fill( avalanches_k.begin(), avalanches_k.end(), 0.0 );
	fill( k_vs_t.begin(), k_vs_t.end(), 0.0 );

	for( size_t smpl = 0; smpl < num_samples; smpl++ ) {
	
		cout << "sample " << smpl+1 << endl;
		Lat1->activateSites();

		t = 0;
		k = 0;
		k_last = 0;

		/* Evacuate the lattice, recording the history avalanche sizes
		   in the vector avalanches. */

		while( Lat1->getNumActive() > 0 ) {
			removeActiveSite( Lat1, CullingQ, t, k, rng );
			a_size = cullSites(Lat1, CullingQ, mcull) - 1;

			k_vs_t[t] += k;
			avalanches_t[t] += a_size;
			for( size_t j = k_last+1; j < k; j++ ) {
				avalanches_k[j] += -1;
			}
			avalanches_k[k] += a_size;
			k_last = k;

		} /* while */

		if( k > k_max ) k_max = k;
		if( t > t_max ) t_max = t;
		
		t_e += t;
		k_e += k;

	} /* smpl */

	for( size_t i = 0; i < num_sites+1; i++ ) {
		avalanches_t[i] /= num_samples;
		avalanches_k[i] /= num_samples;
		k_vs_t[i] /= num_samples;
	}

	t_e /= num_samples;
	k_e /= num_samples;

	cout << "\nPrinting data to files.";

	stringstream filename;
	filename << "avalanche_kt_" << mcull;
	gi.printType(filename); /* Append graph info to file name */
	filename << ".dat";
	ofstream avt_file( filename.str().c_str(), ios::app );
	avt_file.setf(ios::left, ios::adjustfield);
	avt_file << "# ";
    avt_file << filename.str() << endl;
	
	/* print avalanche vs t, k, k_vs_t */
    gi.printGraphInfo(avt_file, num_sites, mcull, num_samples);
	avt_file << "# t_e = " << t_e/num_sites << " k_e = " << k_e/num_sites << endl;
	avt_file <<  "# t or k\t" <<  "avalanche_t\t" << "avalanche_k\t" << "k_vs_t\t" << endl;
	for( size_t i = 1; i <= (size_t)k_e; i++ ) {
		avt_file << (double)i/num_sites << '\t' << avalanches_t[i];
		avt_file << '\t' << avalanches_k[i] << '\t' << k_vs_t[i]/num_sites << endl;
	} 
	avt_file << endl;

	cout << endl;

	return 0;
}
