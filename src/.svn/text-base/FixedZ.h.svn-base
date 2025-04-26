/* FixedZ.h
 * This class is meant to be a superclass random graphs.  It is meant for site
 * percolation only. A different scheme for keeping track of bonds is necessary
 * for bond-percolation.
 *
 * The neighbor retrieval scheme is a neighbor table.
 *
 * 08/01/05 Created this as a subclass of Lattice.
*/

#ifndef FIXEDZ_H
#define FIXEDZ_H

#include "LatticeRandom.h"

class FixedZ: public LatticeRandom
{
    public:
         FixedZ(size_t len, size_t id = 0, size_t z = 6) : LatticeRandom( len, id )
         {
             Setup(z);
             generateBonds();
         }

         virtual ~FixedZ() 
         {
             Cleanup();
         }

         virtual void generateBonds();
         virtual void activateSites();

         virtual size_t getNumNeighbors( size_t); 

    private:
         void Setup(size_t z);
         void Cleanup();
         size_t *coordination;
         size_t *choose_i;
};

#endif
