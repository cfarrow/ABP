#include "LatticeRandom.h"
#include <queue>

size_t LatticeRandom::getNumActiveNeighbors(size_t i) {
    size_t count = 0;
	for(size_t j : getNbrs(i))
    {
        if( isActive( j ) ) count++;
    }
    return count;
}

/* The bond scheme for random, or more appropriately, non-regular graphs is to
 * have a neighbor table. This is memory inefficient, since each bond exists
 * twice in the table, e.g., i is a neighbor of j, which is a neighbor of i.
 * However, it is fast and simple This scheme does not work for removing bonds.
 * It only works for site-percolation.
 */
void LatticeRandom::generateBonds() {
    //IMPLEMENT IN DERIVED CLASS
}

void LatticeRandom::labelClusters(size_t a) {

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

			for(size_t s2 : getNbrs(s1)) {
				if( isActive(s2) == a && cluster_label[s2] == num_sites ) {
									 
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