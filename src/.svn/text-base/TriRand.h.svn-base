/* TriRand.cpp 02/13/07 - initial version of file
 *
 * A 'TriRand' lattice is a triangular lattice where a fraction, f, of bonds are
 * removed and reconnected at long range. Note that periodic boundary conditions
 * are enforced in all directions, as this is essentially a random graph.
*/

#ifndef TRIRAND_H
#define TRIRAND_H

#include "LatticeRandom.h"

class TriRand: public LatticeRandom
{
    public:
        TriRand(size_t len, size_t id = 0, double fb = 0) : LatticeRandom( len*len, id )
        {
            Setup(fb);
            generateBonds();
        }
        
        virtual ~TriRand() 
        {
            Cleanup();
        }
        
        virtual void generateBonds();
        virtual void activateSites();
        
        virtual size_t getNumNeighbors(size_t) {
            return max_neighbors;
        }

    private:
        double f;
        size_t nb; // # of bonds to be swapped out
        int *bondlist;
        void Setup(double fb);
        void Cleanup();
};

#endif
