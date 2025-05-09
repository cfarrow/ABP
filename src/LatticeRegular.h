/* LatticeRegular.h
 * This class is meant to be a superclass for all regular graphs.  It's neighbor
 * retrieval scheme does not work with random graphs.
 *
 * 03/14/07 Created this as a super-class for regular lattices.
*/

#ifndef LATTICEREGULAR_H
#define LATTICEREGULAR_H

#include "Lattice.h"
#include "Neighbors.h"


class LatticeRegular: public Lattice
{
    public:
        LatticeRegular(size_t size, size_t id=0, size_t nn=0) 
        : Lattice( size, id ), num_neighbors(nn) {
            nbrs = new size_t[num_neighbors];
            last = num_sites;
        }

        virtual ~LatticeRegular(){
            delete[] nbrs;
        }

        /* Get an object that iterates over neighbors of a given site. */
        Neighbors getNbrs(size_t);
        
        /* virtual functions */
        virtual void activateSites() {
            Lattice::activateSites();
            last = num_sites;
        }
        virtual size_t getNumNeighbors(size_t) = 0; 
        virtual size_t getNbr(size_t, size_t);
        virtual size_t getNumActiveNeighbors(size_t);

    protected:
         size_t num_neighbors;
         size_t last;
         size_t* nbrs;
         virtual size_t setNbrs(size_t) = 0;
};

#endif
