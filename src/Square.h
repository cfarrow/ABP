/* Square.h
 * 08/01/05 Created this as a subclass of LatticeRegular2d.
*/

#ifndef SQUARE_H
#define SQUARE_H

#include "LatticeRegular2d.h"


template<short pbc>
class Square: public LatticeRegular2d
{
    public:
        Square(size_t len, size_t id = 0) : LatticeRegular2d(len, id, 4) {}
        virtual ~Square() {}
        virtual size_t getNumNeighbors(size_t);

    protected:
        virtual void setNbrs(size_t);
};

#endif
