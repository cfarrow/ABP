/* UJack.h
 * 08/01/05 Created this as a subclass of LatticeRegular2d.
*/

#ifndef UJACK_H
#define UJACK_H

#include "LatticeRegular2d.h"


template<short pbc>
class UJack: public LatticeRegular2d
{
    public:
        UJack(size_t len, size_t id = 0) : LatticeRegular2d(len, id, 8) {}
        virtual ~UJack() {}
        virtual size_t getNumNeighbors(size_t);

    protected:
        virtual size_t setNbrs(size_t);

};
#endif
