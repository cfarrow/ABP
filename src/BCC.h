/* BCC.h
 * 03/06/07 Created this as a subclass of LatticeRegular3d.
*/

#ifndef BCC_H
#define BCC_H

#include "LatticeRegular3d.h"

template<short pbc>
class BCC: public LatticeRegular3d
{
    public:
        BCC(size_t len, size_t id=0) : LatticeRegular3d(len, id, 8) {}
        virtual ~BCC() {}
        virtual size_t getNumNeighbors(size_t i);

    protected:
        virtual void setNbrs(size_t);
};
#endif
