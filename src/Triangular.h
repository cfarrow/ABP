/* Triangular.h
 * 08/01/05 Created this as a subclass of LatticeRegular2d.
*/

#ifndef TRIANGULAR_H
#define TRIANGULAR_H

#include "LatticeRegular2d.h"


template<short pbc>
class Triangular: public LatticeRegular2d
{
    public:
        Triangular(size_t len, size_t id = 0) : LatticeRegular2d(len, id, 6) {}
        virtual ~Triangular() {}
        virtual size_t getNumNeighbors(size_t);

    protected:
        virtual void setNbrs(size_t);
};

#endif
