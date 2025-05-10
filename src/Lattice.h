/* Lattice.h
   02/03/04 An attempt at creating some classes for lattices.
   02/05/04 Modified Member functions
   04/19/04 Messed with it some more. Needed to make site information more
        accessible.  Decided to leave most functions out of the class
        constructor. This makes it general so that optional data can be used
        only if neccessary.
   04/27/04 This "binary" version is for very large lattices. Each site of the
        lattice is a bit, either '0' (inactive) or '1' (active). This should
        increase storage capacity by a factor of 32, which will allow testing of
        larger lattices. Unfortunately, the method of checking the spanning
        cluster still needs an array of cluster labels. This is the largest data
        structure that needs to be accounted for.
   08/11/04 Reworked specifically for ABP. Periodic BC enforced on all sides. No
       BCC due to difficulty.
   10/06/04 Added cluster labeling.
   11/18/04 Added spanning.
   08/01/05 Modified class so that specific lattice types can be subclassed.
       Note that the Lattice class still takes care of most of the dirty work.
       The Lattice constructor now takes the number of sites and an id. It used
       to take the lattice length and the lattice type. The lattice id can be
       used to indicate the latice type. See ABPprocs.cpp and ABPprocs.h for how
       this is done.
*/

#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <algorithm>
#include "Neighbors.h"


/* ============================= Utilities ================================== */

// Size of a "word" representing multiple sites, one bit per site.
const size_t word_size = 8 * sizeof (size_t);

// Get a bit mask for manipulating the state of a location within a word
size_t bit_mask(size_t shift);

// Get the bit value at a specified site within an array of words
bool get_site_value(size_t site_index, size_t* words);

// Set the bit value at a specified site within an array of words
void set_site_value(size_t site_index, bool value, size_t* words);

/* ============================ Lattice Class Functions ===================== */

class Lattice 
{
     public:
          /* The second entry is for an optional type identifier, specified in
          * the program. 
          */
          Lattice(size_t nsites = 32, size_t id = 0);
          virtual ~Lattice();

          /* Functions that act on whole lattice */
          virtual void activateSites(); 
          virtual void labelClusters(size_t=1);
          size_t getLength() {return length;}   	 
          size_t getNumSites() {return num_sites;}
          size_t getNumActive() {return num_active_sites;}
          size_t getNumPresent() {return num_present_sites;}
          size_t getNumClusters() {return Clusters.size();}
          size_t getMaxMass() {return (Clusters.empty() ? 0 : *(max_element(Clusters.begin(), Clusters.end()))); }

          /* Functions that act on single sites */
          bool isActive(size_t);
          bool isPresent(size_t);
          /* the following do not do any checking. It is up to the algorithm. 'if' statements slow things down! */
          void setActiveLevel(size_t, bool);
          void setPresentLevel(size_t, bool);
          /* get cluster label of site */
          size_t getClusterLabel(size_t i) {return cluster_label[i]; }
          /* get size of cluster, as referenced by its label */
          size_t getClusterSize(size_t i) {return Clusters[i]; }
          size_t getNumActiveNeighbors(size_t);
          size_t getDims() {return dims;}

          /* ---------------------------------------------------------------*\
          Check to see if the largest cluster is a spanning cluster.
          Looks to see if a site from both the top and bottom row (plane)
          are on the largest cluster. This must be defined in the specific
          lattice-type class (see Triangular.cpp for an example.)
          \* ---------------------------------------------------------------*/
          virtual bool isSpanning(size_t = 1) = 0;

          /* virtual function - see notes below */
          virtual size_t getNumNeighbors(size_t) = 0; 
          virtual size_t getNbr(size_t, size_t) = 0;
          virtual Neighbors getNbrs(size_t, bool) = 0;
          virtual Neighbors getNbrs(size_t) = 0;

     protected: /* These are direcly callable/mutable by derived classes */
          size_t length;			/* length of lattice - if applicable */
          size_t dims;             /* the dimensionality of the lattice */
          size_t graph_id;         /* number id of graph, could represent type */
          size_t num_sites;		/* number of sites	*/
          size_t num_words;
          size_t num_active_sites;	/* number of active sites */
          size_t num_present_sites;	/* number of present sites */

          size_t *active;				/* active indicator */
          size_t *present;			/* present indicator */
          size_t *cluster_label; /* array of cluster labels. Initialized in labelClusters */
          std::vector < size_t > Clusters;
};
#endif