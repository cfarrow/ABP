/* Lattice.cpp
   02/03/04 Member functions for Lattice.h
   02/05/04 Modified member functions.
   04/19/04 Reworked for speed.
   04/27/04 This "binary" version is for very large lattices. Each site of the
	    lattice is a bit, either '0' (inactive) or '1' (active). This should
	    increase storage capacity by a factor of 32, which will allow testing
	    of larger lattices. 
   08/11/04 Reworked specifically for ABP. Periodic BC enforced on all sides.
   		No BCC yet due to difficulty in coding neighbor function.
   10/06/04 Added cluster labeling.
   11/18/04 Added spanning cluster check.

*/

#include <iostream>
#include <queue>

#include "Lattice.h"


/* ============================= Utilities ================================== */

inline size_t bit_mask(size_t shift) {
	// mask with 'word_size' bits, all zero except the bit at location 'shift'
	return size_t{1} << shift;
}

inline bool get_site_value(size_t site_index, size_t* words) {
	size_t shift = site_index % word_size;
	size_t word_num = (site_index - shift) / word_size;
	size_t word = words[word_num];

	return word & bit_mask(shift);
}

inline void set_site_value(size_t site_index, bool value, size_t* words) {
	size_t shift = site_index % word_size;
	size_t word_num = (site_index - shift) / word_size;

	if(value) {
		words[word_num] |= bit_mask(shift);
	}
	else {
		words[word_num] &= ~bit_mask(shift);
	}
}

/* ============================ Lattice Class Functions ===================== */

Lattice::Lattice(size_t nsites, size_t id) {
		 
         graph_id = id;
		 num_active_sites = 0;
		 num_present_sites = 0;
		 cluster_label = 0;

         /* Make sure the number of sites is divisible by word_size */
         num_sites =  ( nsites / word_size ) * word_size;

         /* num_sites must be a multiple of word_size!*/
		 num_words = num_sites / word_size;
		 /* allocate memory for essential class data */
		 active = new size_t [num_words];
		 present = new size_t [num_words];
		 if( active == NULL || present == NULL ) {
             std::cerr << "Memory allocation error!\n";
			 exit(1);
		 }

		 return;
}

bool Lattice::isActive(size_t i) {
	return get_site_value(i, active);
}

bool Lattice::isPresent(size_t i) {
	return get_site_value(i, present);
}

void Lattice::setActiveLevel(size_t i, bool new_level) {

	bool old_level = get_site_value(i, active);
	set_site_value(i, new_level, active);

	/* Update the active neighbors.  */
	if( new_level && !old_level ) { 
		num_active_sites++;
	}
	else if( !new_level && old_level ) { 
		num_active_sites--; 
	}

	return;
}

void Lattice::setPresentLevel(size_t i, bool new_level) {

	bool old_level = get_site_value(i, present);
	set_site_value(i, new_level, present);

	/* Update the active neighbors.  */
	if( new_level && !old_level ) { 
		num_present_sites++;
	}
	else if( !new_level && old_level ) { 
		num_present_sites--; 
	}

	return;
}


Lattice::~Lattice() {
	if( cluster_label != 0 ) delete [] cluster_label;
	delete [] active;
	delete [] present;
}

/* Set all sites to active and present */
/* May require rewriting in a derived class
 * if such class needs to know if the lattice
 * has just been activated.
 */
void Lattice::activateSites(void) {
	 
	num_words = num_sites / word_size;
	static size_t ones = ~size_t{0}; // Every bit is 1
	for( size_t i = 0; i < num_words; i++ ) {
			active[i] = ones;
			present[i] = ones;
	}
	num_active_sites = num_sites;
	num_present_sites = num_sites;
}

/* ---------------------------------------------------------------*\

   Find the clusters in the lattice and label the sites
   according to their cluster. Must relax neighboring conditions
   since periodic BC may be used. This comes down
   to checking if a site is in x = 0 or x = length -1 plane.
   If it is, then we ignore it's neighbor if it is
   in the periodically neighboring row (plane). For the 2D lattices
   we check y-values. In 3D lattices, we check z-values as well.
   For random graphs, this condition is left out.

\* ---------------------------------------------------------------*/
void Lattice::labelClusters(size_t a) {

	/* initialize the cluster_label array if it has not been already */
	if( cluster_label == 0 ) {
		cluster_label = new size_t[num_sites];
	}

	/* a indicates which sites to label. a = 1 labels active sites.
	 * a = 0 labels inactive sites. a > 1 is not allowed.  */
	if( a > 1 ) { a = 1; }

    std::queue< size_t > CQ; /* The cluster queue */
	size_t s1; /* for holding sites */
	size_t num_clusters = 0;

	/* initialize to large label */
	for( size_t i=0; i < num_sites; i++ ) {
		cluster_label[i] = num_sites;
	}

	Clusters.erase( Clusters.begin(), Clusters.end() );
  
	for( size_t i0=0; i0 < num_sites; i0++ ){
		/* if active != a or already labeled, continue */
		if( isActive(i0) != a || cluster_label[i0] != num_sites ) continue;
    
		/* start a new cluster */
		Clusters.push_back( 0 );
    
		/* label this site and count its mass */
		cluster_label[i0] = num_clusters;
		Clusters[num_clusters]++;
    
		/* send to queue */
		CQ.push( i0 );
    
		while ( ! CQ.empty() ) { /* while queue is nonempty */
      
			/* extract site from queue */ 
			s1 = CQ.front();
			CQ.pop();

			/* Whenever we check a site, find its row (plane). If it's row is
			   0, it is in the top row. If its row is length-1, then it is in
			   the bottom row. If the difference in the planes of the
			   two sites is length-1, then they must be in in the top and 
			   bottom planes, respectivley.
			 */
			for(size_t s2 : getNbrs(s1)) {
				if( isActive(s2) == a && cluster_label[s2] == num_sites && 
					(size_t) abs(static_cast<int>(s1%length) - static_cast<int>(s2%length)) != length-1 ) {
									 
					/* label this site and count its mass */
					cluster_label[s2] = num_clusters;
					Clusters[num_clusters]++;
									 
					/* send to queue */
					CQ.push( s2 );

				} /* if */

			} /* for j */
		
		} /* while */

    num_clusters++; /* increase cluster label since all sites on a given cluster were visited */

  } /* for i0 */

  return;
}