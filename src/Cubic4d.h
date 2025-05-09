/* Cubic4d.h
 * 03/15/07 Created this as a subclass of LatticeRegular4d.
*/

#ifndef CUBIC4D_H
#define CUBIC4D_H

#include "LatticeRegular4d.h"

template<short pbc>
class Cubic4d: public LatticeRegular4d
{
    public:
        Cubic4d(size_t len, size_t id=0) : LatticeRegular4d(len, id, 8)  {}
        virtual ~Cubic4d() {}
        virtual size_t getNumNeighbors(size_t i);

    protected:
        virtual size_t setNbrs(size_t);
};
#endif
