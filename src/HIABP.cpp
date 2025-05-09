/* -------------------------------------------------------------------------------*\ 
 04/16/04 Program takes a lattice and finds the EXACT BP transition point using 
 the avalanching algorithm and a half interval search. The procedure starts by
 randomly labeling the sites. It then completes a full avalanche of the lattice
 and tags each site with the avalanche number it got caught up in. It then picks
 two reasonable bounds and does a half interval search to find the exact avalanche
 cluster that disconnected the lattice by replacing sites via their avalanche
 cluster label. Sites that are being moved back and forth between instances
 are kept in two stacks for easy retrieval.
 05/03/05 Rewrote to be compatible with new class hierarchy for Lattice.
\* -------------------------------------------------------------------------------*/

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <queue>
#include <sstream>
#include <stack>

#include "Lattice.h"  /* compile with "Lattice.cpp" */
#include "ABPprocs.h" /* compile with "ABPprocs.cpp" */
#include "rand.h"
using namespace std;


/* Prototypes */
int cullSitesII(Lattice *L, queue< size_t > &Q, size_t m_val,
		stack< size_t > &RemS, size_t Avalanche[], size_t &a_label);
void restoreSites( Lattice *L, stack< size_t> &RemS, 
		stack< size_t> &ResS, size_t Avalanche[], size_t &t_lo );
void removeSites( Lattice *L, stack< size_t> &RemS, 
		stack< size_t> &ResS, size_t Avalanche[], size_t &t_hi );
void exit_and_usage(void);

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
	stack< size_t, deque<size_t> > rem_stack; /* stack of removed sites */
	stack< size_t, deque<size_t> > res_stack; /* stack of restored sites */
	
    // Mon Jan  9 2006 - CLF
    // Changed the 'var' variables to 'unc' variables. They give uncertainty in
    // the mean of the measured value.
	double tc_avg = 0.0; 
	double tc_unc = 0.0; 
	double tf_avg = 0.0; 
	double tf_unc = 0.0; 
	double pc_avg = 0.0;
	double pc_unc = 0.0;
	double P_avg = 0.0;
	double P_unc = 0.0;

	bool is_spanning;
	size_t k = 0; /* label of the next site to be removed (or just removed) */
	size_t t = 0;
	vector< size_t > k_vs_t(num_sites);
	vector< double > c_array(num_samples);
	vector< double > tc_array(num_samples);
	vector< double > tf_array(num_samples);
	vector< double > P_array(num_samples);

	size_t *avalanche_label = new size_t [num_sites];
	if( avalanche_label == NULL ) {
		cerr << "Memory allocation error!" << endl;
		exit(1);
	}

	bounded_rng_type rng = makeRNG(7, 0, num_sites-1);

	for( size_t smpl = 0; smpl < num_samples; smpl++ ) {
		k = 0;
		t = 0;
	
		Lat1->activateSites();	
		cout << "sample " << (smpl+1) << endl;

		while(! res_stack.empty() ) res_stack.pop();
		while(! rem_stack.empty() ) rem_stack.pop();

		/* Evacuate the lattice, recording the history of removed sites
		   in rem_stack and in avalanche_label. avalanche_label is the
		   t-value at the time of the avalanche. 
		 */
		   
		while( Lat1->getNumActive() > 0 ) {
			removeActiveSite(Lat1, CullingQ, t, k, rng);
			cullSitesII(Lat1, CullingQ, mcull, rem_stack, avalanche_label, t); 
			k_vs_t[t] = k;

		} /* while */
        tf_array[ smpl ] = (double)t/num_sites;

		/* The number of clusters is recorded and each site now has a
		   culling cluster label attached to it. Now take some reasonable
		   bounds for t, the number of steps to disconnect the lattice and
		   do a half interval search between the bounds for the number of
		   the culling cluster that disconnected the lattice. Note that
		   the higher values of t are near the top of the stack. num_steps
		   is the number of half-interval steps used to locate t. 
		 */
		
		size_t t_high = t;
		size_t t_low = 1;
		size_t t_estimate, t_temp;
		size_t num_steps;

		t_estimate = t_low;
		restoreSites( Lat1, rem_stack, res_stack, avalanche_label, t_estimate );
		Lat1->labelClusters();
		is_spanning = Lat1->isSpanning();

		/* The number of half-interval steps */
		num_steps = (size_t) ( ceil( log((double) num_sites) / log(2.0) ) - 1 );

		/* Note that t_estimate is the the smallest culling cluster
		   label present in res_stack. All sites in res_stack
		   are active. When the for loop terminates, if the
		   stable cluster is spanning, then it is the last spanning
		   cluster and thus t_estimate = t. In addition, the
		   top site in rem_stack is the site (that maps to) k+1 since
		   the removal of this site disconnects the lattice.
		   If it is not spanning, then it is the first broken cluster.
		 */
		for( size_t i = 0; i < num_steps; i++ ) {
			
			/* if the cluster is spanning (t is too low, i.e. not enough 
			   sites have been removed) raise t_estimate and
			   raise t_low, then remove some sites. */
			if(is_spanning) {
				t_temp = t_estimate;
				t_estimate = (size_t)( ceil( 0.5*(t_high+t_estimate) ) );
				t_low = t_temp;
				removeSites( Lat1, rem_stack, res_stack, avalanche_label, t_estimate );
			}

			/* if the cluster is not spanning (t is too high, i.e. too many 
			   sites have been removed) lower t_estimate and
			   lower t_high, then add some sites. */
			else {
				t_temp = t_estimate;
				t_estimate = (size_t) floor( 0.5*(t_estimate+t_low) ) ;
				t_high = t_temp;
				restoreSites( Lat1, rem_stack, res_stack, avalanche_label, t_estimate );
			}

			Lat1->labelClusters();
			is_spanning = Lat1->isSpanning();	
		}/* find last spanning cluster */
		
		while(!is_spanning) {
			t_estimate--;
			restoreSites( Lat1, rem_stack, res_stack, avalanche_label, t_estimate );
			Lat1->labelClusters();
			is_spanning = Lat1->isSpanning();
		}

		/* res_stack now holds the last spanning configuration of the lattice. 
		   The top value of the stack is the last site that, when removed,
		   will lead to the collapse of the giant cluster. 
		 */
		
        tc_array[ smpl ]= (double) t_estimate/num_sites;
		P_array[ smpl ] = (double) Lat1->getMaxMass()/num_sites;
		c_array[ smpl ] = 1 - (double) k_vs_t[t_estimate]/num_sites;
		
	} /* smpl */

	P_avg = accumulate( P_array.begin(), P_array.end(), 0.0)/num_samples;
	pc_avg = accumulate( c_array.begin(), c_array.end(), 0.0 )/num_samples;
	tc_avg = accumulate( tc_array.begin(), tc_array.end(), 0.0 )/num_samples;
	tf_avg = accumulate( tf_array.begin(), tf_array.end(), 0.0 )/num_samples;

	pc_unc = 0;
	P_unc = 0;
    tc_unc = 0;
    tf_unc = 0;
	for( size_t i = 0; i < num_samples; i++) {
			 pc_unc += pow( pc_avg - c_array[i], 2.0);
			 P_unc += pow( P_avg - P_array[i], 2.0);
			 tc_unc += pow( tc_avg - tc_array[i], 2.0);
			 tf_unc += pow( tf_avg - tf_array[i], 2.0);
	}
	if(num_samples > 1) {
		pc_unc = sqrt( pc_unc/(num_samples - 1) );
		tc_unc = sqrt( tc_unc/(num_samples - 1) );
		tf_unc = sqrt( tf_unc/(num_samples - 1) );
		P_unc = sqrt( P_unc/(num_samples - 1) );
	}
    /* Calculate the uncertainty in the mean. */
    pc_unc /= sqrt(num_samples);
    tc_unc /= sqrt(num_samples);
    tf_unc /= sqrt(num_samples);
    P_unc  /= sqrt(num_samples);

	/* print the outcome to screen */
	cout.setf(ios::left, ios::adjustfield);
	cout << setw(8) << "#L" << setw(12) << "samples"; 
    cout << setw(12) << "P" << setw(12) << "P_unc";
	cout << setw(12) << "pc" << setw(12) << "pc_unc"; 
	cout << setw(12) << "kc" << setw(12) << "kc_unc"; 
	cout << setw(12) << "tc" << setw(12) << "tc_unc"; 
	cout << setw(12) << "tf" << setw(12) << "tf_unc" << endl;
	cout << setw(8) << Lat1->getLength() << setw(12) << num_samples;
    cout << setw(12) << P_avg << setw(12) << P_unc;
	cout << setw(12) << pc_avg << setw(12) << pc_unc;
	cout << setw(12) << 1 - pc_avg << setw(12) << pc_unc;
	cout << setw(12) << tc_avg << setw(12) << tc_unc;
	cout << setw(12) << tf_avg << setw(12) << tf_unc << endl;

	/* print the outcome to file */
	stringstream filename;
	filename << "HIABP_" << mcull;
	gi.printType(filename); /* Append graph info to file name */
    filename << ".dat";
	bool new_file = false;
	ifstream temp_file( filename.str().c_str(), ios::in );
	if ( !temp_file ) new_file = true;
	temp_file.close();
	ofstream out_file( filename.str().c_str(), ios::app );
	out_file.setf(ios::left, ios::adjustfield);
	if( new_file ) {
        out_file << setw(8) << "#L" << setw(12) << "samples" ;
        out_file << setw(12) << "P" << setw(12) << "P_unc";
        out_file << setw(12) << "pc" << setw(12) << "pc_unc"; 
        out_file << setw(12) << "kc" << setw(12) << "kc_unc"; 
        out_file << setw(12) << "tc" << setw(12) << "tc_unc"; 
        out_file << setw(12) << "tf" << setw(12) << "tf_unc" << endl;
	}
	out_file << setw(8) << Lat1->getLength() << setw(12) << num_samples;
    out_file << setw(12) << P_avg << setw(12) << P_unc;
	out_file << setw(12) << pc_avg << setw(12) << pc_unc;
	out_file << setw(12) << 1 - pc_avg << setw(12) << pc_unc;
	out_file << setw(12) << tc_avg << setw(12) << tc_unc;
	out_file << setw(12) << tf_avg << setw(12) << tf_unc << endl;

	return 0;
}

/* This culls sites just as in the ABPprocs.cpp file, but it adds removed sites to
   the rem_stack and also labels the sites according to what culling_cluster they
   are a part of.  */
int cullSitesII(Lattice *L, queue< size_t > &Q, size_t m_val, 
	stack< size_t > &RemS, size_t Avalanche[], size_t &a_label) {

	size_t s1, n_cul = 0;
  
	while (!Q.empty()) {		/* while queue still has sites  */
		
		s1 = Q.front();
	    Q.pop();
		RemS.push(s1);		/* put the site on the removed stack */
		Avalanche[s1] = a_label;	/* label the site with its culling cluster */
		L->setActiveLevel(s1,0);
		n_cul++;
      
		/* Now check all neighbor sites of the recently culled site and send them to the
		   culling queue if they don't meet the culling condition.                         */
	    for(size_t s2 : L->getNbrs(s1)) {
			if ( L->isActive(s2) && L->getNumActiveNeighbors(s2) < m_val ) {
				Q.push(s2);
				L->setActiveLevel(s2,0);
			}/* if */
		}/* for j */
	}/* while */

	return n_cul;
}

/* Restore the sites with culling avalanche label up through t_lo to active */
void restoreSites( Lattice *L, stack< size_t> &RemS, 
		stack< size_t> &ResS, size_t Avalanche[], size_t &t_lo ) {

	size_t s1;

	do {
		s1 = RemS.top();
		RemS.pop();
		ResS.push(s1);
		L->setActiveLevel(s1,1);
		if( RemS.size() == 0 ) break;

	} while( Avalanche[ RemS.top() ] >= t_lo );
}

/* Make inactive the sites with avalanche cluster label up to (but not including) t_hi */
void removeSites( Lattice *L, stack< size_t> &RemS, 
		stack< size_t> &ResS, size_t Avalanche[], size_t &t_hi ) {

	size_t s1;

	do {
		s1 = ResS.top();
		ResS.pop();
		RemS.push(s1);
		L->setActiveLevel(s1,0);
		if( ResS.size() == 0 ) break;
		
	} while( Avalanche[ RemS.top() ] < t_hi );
}

void exit_and_usage(void)  {

	cout << "usage: " << "HIABP" << " <type> <length> <m> <samples>\n";
	exit(1);
}

