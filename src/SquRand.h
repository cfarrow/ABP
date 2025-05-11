/* SquRand.cpp 04/13/07 - initial version of file
 *
 * A 'SquRand' lattice is a square lattice where a fraction, f, of bonds are
 * removed and reconnected at long range. Note that periodic boundary conditions
 * are enforced in all directions, as this is essentially a random graph.
*/

#ifndef SQURAND_H
#define SQURAND_H

#include "LatticeRandom.h"

class SquRand: public LatticeRandom
{
    public:
         SquRand(size_t len, size_t id = 0, double fb = 0) : LatticeRandom( len*len, id )
         {
             Setup(fb);
             generateBonds();
         }

         virtual ~SquRand() 
         {
             Cleanup();
         }

         virtual void generateBonds();
         virtual void activateSites();

         virtual size_t getNumNeighbors(size_t) {return max_neighbors;}
         size_t getDims() {return 2;}

    private:
         double f;
         size_t nb; // # of bonds to be swapped out
         int *bondlist;
         void Setup(double fb);
         void Cleanup();
};

#endif
