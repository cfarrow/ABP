/* -------------------------------------------------------------------------------*\ 
 03/15/07 Record the maximum avalalanche.
\* -------------------------------------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <cmath>
#include <queue>
#include "Lattice.h"  /* compile with "Lattice.cpp" */
#include "ABPprocs.h" /* compile with "ABPprocs.cpp" */
#include "rand.h"
using namespace std;

int main(int argc, char **argv) {

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

	queue< size_t > CullingQ; /* The culling queue. */
	
	double a_max_avg = 0.0; 
	double a_max_unc = 0.0; 
	double t_max_avg = 0.0; 
	double t_max_unc = 0.0; 
	double k_max_avg = 0.0; 
	double k_max_unc = 0.0; 
    double t_e = 0.0;
    double k_e = 0.0;

	bool is_spanning;
    size_t a_size;
	size_t k = 0; /* label of the next site to be removed (or just removed) */
	size_t t = 0;

	size_t *a_max = new size_t [num_samples];
	double *t_max = new double [num_samples];
	double *k_max = new double [num_samples];
    for(size_t i=0; i<num_samples; ++i)
    {
        a_max[i] = 0;
        t_max[i] = 0.0;
        k_max[i] = 0.0;
    }

	bounded_rng_type rng = makeRNG(7, 0, num_sites-1);
	for( size_t smpl = 0; smpl < num_samples; smpl++ ) {
	
		cout << "sample " << smpl+1 << endl;
		Lat1->activateSites();

		t = 0;
		k = 0;

		while( Lat1->getNumActive() > 0 ) {
			removeActiveSite(Lat1, CullingQ, t, k, rng);
			a_size = cullSites(Lat1, CullingQ, mcull) - 1;
            if(a_size > a_max[smpl]) {
                a_max[smpl] = a_size;
                t_max[smpl] = t;
                k_max[smpl] = k;
            }
		} /* while */

		t_e += t;
		k_e += k;

	} /* smpl */

    // Calculate averages
	t_e /= num_samples*num_sites;
	k_e /= num_samples*num_sites;

    for(size_t i=0; i<num_samples; ++i) {
        a_max_avg += a_max[i];
        // Convert (t,k) to (tau,kappa)
        t_max[i] /= num_sites;
        t_max_avg += t_max[i];
        k_max[i] /= num_sites;
        k_max_avg += k_max[i];
    }
    a_max_avg /= num_samples;
    t_max_avg /= num_samples;
    k_max_avg /= num_samples;

    // Calculate uncertainties
	for( size_t i = 0; i < num_samples; i++) {
         a_max_unc += pow( a_max_avg - a_max[i], 2.0);
         t_max_unc += pow( t_max_avg - t_max[i], 2.0);
         k_max_unc += pow( k_max_avg - k_max[i], 2.0);
	}
	if(num_samples > 1) {
        a_max_unc = sqrt( a_max_unc/(num_samples -1) );
        t_max_unc = sqrt( t_max_unc/(num_samples -1) );
        k_max_unc = sqrt( k_max_unc/(num_samples -1) );
	}

    /* Calculate the uncertainty in the mean. */
    a_max_unc /= sqrt(num_samples);
    t_max_unc /= sqrt(num_samples);
    k_max_unc /= sqrt(num_samples);

	/* print the outcome to screen */
	cout.setf(ios::left, ios::adjustfield);
	cout << setw(8) << "#L" << setw(12) << "samples"; 
	cout << setw(12) << "a_max_avg" << setw(12) << "a_max_unc"; 
	cout << setw(12) << "t_max_avg" << setw(12) << "t_max_unc"; 
	cout << setw(12) << "k_max_avg" << setw(12) << "k_max_unc"; 
	cout << setw(12) << "t_e" << setw(12) << "k_e" << endl;
	cout << setw(8) << Lat1->getLength() << setw(12) << num_samples;
	cout << setw(12) << a_max_avg << setw(12) << a_max_unc;
	cout << setw(12) << t_max_avg << setw(12) << t_max_unc;
	cout << setw(12) << k_max_avg << setw(12) << k_max_unc;
	cout << setw(12) << t_e << setw(12) << k_e << endl;

	/* print the outcome to file */
	stringstream filename;
	filename << "Avalanche_max_" << mcull;
	gi.printType(filename); /* Append graph info to file name */
    filename << ".dat";
	bool new_file = false;
	ifstream temp_file( filename.str().c_str(), ios::in );
	if ( !temp_file ) new_file = true;
	temp_file.close();
	ofstream out_file( filename.str().c_str(), ios::app );
	out_file.setf(ios::left, ios::adjustfield);
	if( new_file ) {
        out_file << setw(8) << "#L" << setw(12) << "samples"; 
        out_file << setw(12) << "a_max_avg" << setw(12) << "a_max_unc"; 
        out_file << setw(12) << "t_max_avg" << setw(12) << "t_max_unc"; 
        out_file << setw(12) << "k_max_avg" << setw(12) << "k_max_unc"; 
        out_file << setw(12) << "t_e" << setw(12) << "k_e" << endl;
	}
	out_file << setw(8) << Lat1->getLength() << setw(12) << num_samples;
	out_file << setw(12) << a_max_avg << setw(12) << a_max_unc;
	out_file << setw(12) << t_max_avg << setw(12) << t_max_unc;
	out_file << setw(12) << k_max_avg << setw(12) << k_max_unc;
	out_file << setw(12) << t_e << setw(12) << k_e << endl;

	return 0;
}
