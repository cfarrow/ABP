/* LatticeRandom.h 
 *
 * This class is a superclass random graphs.  It is meant for site percolation
 * only. A different scheme for keeping track of bonds is necessary for
 * efficient bond-percolation.  The neighbor retrieval scheme is a neighbor
 * table. This should be initialized during the class construction and destroyed
 * at class destruction. This must be implemented in a derived class.  A proper
 * subclass will need to define how to generate and connect bonds. This would
 * take place in a function that gets called every time the lattice is reset.
 *
 * See FixedZ.h for an example.
 *
 * 08/01/05 Created this as a subclass of Lattice.
*/

#ifndef LATTICERANDOM_H
#define LATTICERANDOM_H

#include "Lattice.h"
#include "Neighbors.h"

class LatticeRandom: public Lattice
{
    public:
        LatticeRandom(size_t nsites, size_t id = 0) : Lattice(nsites, id) {}
        virtual ~LatticeRandom() {}
        

        // These functions are specific to this class
        virtual void generateBonds();
        virtual bool isGiant() {
           /* Cluster is giant if it has > 10% of the sites */
           if( getMaxMass() > num_sites/10 ) return true;
           else return false;
        }
        
        Neighbors getNbrs(size_t i) { 
            return Neighbors(neighbor[i], getNumNeighbors(i)); 
        }
        virtual size_t getNbr(size_t i, size_t j) { return neighbor[i][j]; }
        virtual size_t getNumActiveNeighbors(size_t i);
        virtual bool isSpanning(size_t = 1){ return isGiant();}
        virtual void labelClusters(size_t = 1);

    protected:
        size_t max_neighbors; // Maximum # of neighbors for a site.
        size_t *nbr_ref; //Helps make data sequential in memory, will speed up program.
        size_t **neighbor; //The actual neighbor table.
};

#endif
