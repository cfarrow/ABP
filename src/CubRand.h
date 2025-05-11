/* CubRand.cpp 
 *
 * 03/06/07 - initial version of file
 *
 * A 'CubRand' lattice is a cubic lattice where a fraction, f, of bonds are
 * removed and reconnected at long range. Note that periodic boundary
 * conditions are enforced in all directions, as this is essentially a random
 * graph.
*/

#ifndef CUBRAND_H
#define CUBRAND_H

#include "LatticeRandom.h"

class CubRand: public LatticeRandom
{
    public:
        CubRand(size_t len, size_t id = 0, double fb = 0) : LatticeRandom( len*len*len, id )
        {
            Setup(fb);
            generateBonds();
        }
        
        virtual ~CubRand() 
        {
            Cleanup();
        }
        
        virtual void generateBonds();
        virtual void activateSites();
        
        virtual size_t getNumNeighbors(size_t) {return max_neighbors;}
        size_t getDims() {return 3;}

    private:
        double f;
        size_t nb; // # of bonds to be swapped out
        int *bondlist;
        void Setup(double fb);
        void Cleanup();
        size_t getTriNbr(size_t i, size_t nidx);
};

#endif
