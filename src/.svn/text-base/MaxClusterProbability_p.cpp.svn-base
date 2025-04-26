/* -------------------------------------------------------------------------------*\ 
 03/09/07 Records the size of the largest cluster vs k.
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
#include <numeric>
#include "Lattice.h"  /* compile with "Lattice.cpp" */
#include "ABPprocs.h" /* compile with "ABPprocs.cpp" */
#include "utils.h"
#include "MersenneTwister.h"
using namespace std;

/* The global stacks of removed and replaced sites */

/* Prototypes */
static void diluteSites( Lattice *L, queue<size_t> &Q, vector<size_t> &map, size_t &k_min, size_t &k_max );

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
    
    size_t k = 0; /* label of the next site to be removed (or just removed) */
    size_t k_last = 0; /* label of the next site to be removed (or just removed) */

    size_t num_vals = 200; /* The number of k-values to record on */
    double p_int = 1.0/num_vals; /* The p-interval to record on */

    vector< double > P_array_k(num_vals);
    vector< size_t > site_map(num_sites);

    for( size_t i = 0; i < num_vals; ++i )
    {
        P_array_k[i] = 0;
    }

    for( size_t i =0; i < num_sites; ++i )
    {
        site_map[i] = i;
    }

    MTRand(7);
    MTRand mtRNG;
    int num_left, choice;
    size_t s_max;

    for( size_t smpl = 0; smpl < num_samples; smpl++ ) {
        k = 0;
        k_last = 0;
    
        Lat1->activateSites();  
        cout << "sample " << (smpl+1) << endl;

        /* Shuffle the sites in the site_map */
        num_left = num_sites;
        while(num_left > 0)
        {
            choice = static_cast<int>( mtRNG.randExc(num_left) );
            --num_left;
            swap(site_map[choice], site_map[num_left]);
        }

        for(size_t i=0; i < num_vals; ++i) {

            /* Get the current k-label */
            k_last = k;
            k = static_cast<size_t>(i*p_int*num_sites);
            diluteSites(Lat1, CullingQ, site_map, k_last, k);
            cullSites(Lat1, CullingQ, mcull);
            Lat1->labelClusters();
            s_max = Lat1->getMaxMass();
            if(s_max > 0)
            {
                P_array_k[i] += s_max;
            }
            else
            {
                break;
            }
        }
           
    } /* smpl */

    /* Normalize P_array_k */
    for( size_t i = 0; i < num_vals; ++i )
    {
        P_array_k[i]/= num_sites*num_samples;
    }

    /* print the outcome to file */
    stringstream filename;
    filename << "P_max_cluster_p_" << mcull;
	gi.printType(filename); /* Append graph info to file name */
    filename << ".dat";
    ofstream out_file( filename.str().c_str());
    out_file.setf(ios::left, ios::adjustfield);

    gi.printGraphInfo(out_file, num_sites, mcull, num_samples);
    out_file << setw(12) << "# p" << setw(12) << "P_inf" << endl;
    // Record the percolation probability.
    for( int i = num_vals-1; i >= 0; --i )
    {
        out_file << setw(12) << 1.0 - i*p_int;
        out_file << setw(12) << P_array_k[i];
        out_file << endl;
    }

    return 0;
}

static void diluteSites( Lattice* L, queue< size_t > &Q, 
        vector< size_t > &map, size_t &k_min, size_t &k_max ) {

    for( size_t i = k_min; i <= k_max; i++) {
        size_t map_site = map[i];
        Q.push( map_site );
        L->setActiveLevel(map_site, 0);
    }
    return;
}
