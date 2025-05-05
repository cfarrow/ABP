/* Hexagonal.h
 * 01/30/07 Created this as a subclass of LatticeRegular2d.
*/

#ifndef HEXAGONAL_H
#define HEXAGONAL_H

#include "LatticeRegular2d.h"


template<short pbc>
class Hexagonal: public LatticeRegular2d
{
    public:
        Hexagonal(size_t len, size_t id = 0) : LatticeRegular2d(len, id, 3) {}
        virtual ~Hexagonal() {}
        virtual size_t getNumNeighbors(size_t);

    protected:
        virtual void setNbrs(size_t);
};
#endif
