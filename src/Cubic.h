/* Cubic.h
 * 08/01/05 Created this as a subclass of LatticeRegular3d.
*/

#ifndef CUBIC_H
#define CUBIC_H

#include "LatticeRegular3d.h"

template<short pbc>
class Cubic: public LatticeRegular3d
{
    public:
        Cubic(size_t len, size_t id=0) : LatticeRegular3d(len, id, 6) {}
        virtual ~Cubic() {}
        virtual size_t getNumNeighbors(size_t);

    protected:
        virtual void setNbrs(size_t);
};
#endif
