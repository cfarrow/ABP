/* LatticeRegular.h
 * This class is meant to be a superclass for all regular graphs.  It's neighbor
 * retrieval scheme does not work with random graphs.  A proper subclass will
 * only need to define the setNbrs() function.  This means that there is no need
 * for a cpp file.  The derived class may require additional items in the
 * constructor.  See Triangular.h for an example.
 *
 * 03/14/07 Created this as a super-class for regular lattices.
*/

#ifndef LATTICEREGULAR_H
#define LATTICEREGULAR_H

#include "Lattice.h"
#include <vector>

class LatticeRegular: public Lattice
{
    public:
        LatticeRegular(size_t size, size_t id = 0) : Lattice( size, id ) {}
        virtual ~LatticeRegular() {}
        
        /* virtual functions */
        virtual void activateSites() {
            Lattice::activateSites();
            last = num_sites;
        }
        
        /* virtual function - defined below */
        virtual size_t getNumNeighbors( size_t); 
        virtual size_t getNbr(size_t, size_t);
        virtual size_t getNumActiveNeighbors(size_t i);

    protected:
         size_t num_neighbors;
         size_t last;
         std::vector <size_t> nbrs; /* holds the neighbor tags for a given site */
         
         // This function is specific to regular lattices
         void (LatticeRegular::*setNbrs) (size_t); /* populates the nbrs array */

};

typedef void ( LatticeRegular::*LRFptr ) (size_t );

#endif
