/* LatticeRegular4d.h
 *
 * 03/15/07 Created this as a subclass of LatticeRegular.
*/

#ifndef LATTICEREGULAR4D_H
#define LATTICEREGULAR4D_H

#include <cmath>
#include "LatticeRegular.h"


class LatticeRegular4d: public LatticeRegular
{
    public:
        LatticeRegular4d(size_t len, size_t id=0, size_t nn=0) 
        : LatticeRegular(len*len*len*len, id, nn) {
            length = static_cast<size_t>(ceil(pow(num_sites, 0.25)));
            b = length - 1;
        }
        virtual ~LatticeRegular4d() {}

         /* virtual functions */
        virtual bool isSpanning(size_t = 1);

    protected:
        size_t dims{4};
        size_t b; // The boundary of any dimension of the lattice
};
#endif
