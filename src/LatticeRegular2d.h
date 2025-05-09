/* LatticeRegular2d.h
 * This class is meant to be a superclass for 2D regular graphs.  It's neighbor
 * retrieval scheme does not work with random graphs.
 *
 * 08/01/05 Created this as a subclass of Lattice.
 * 03/14/07 Subclassed from LatticeRegular
*/

#ifndef LATTICEREGULAR2D_H
#define LATTICEREGULAR2D_H

#include <cmath>
#include "LatticeRegular.h"


class LatticeRegular2d: public LatticeRegular
{
    public:
        LatticeRegular2d(size_t len, size_t id=0, size_t nn=0)
        : LatticeRegular(len*len, id, nn) {
           length = static_cast< size_t >(ceil(sqrt(num_sites)));
           b = length - 1;
        }
        virtual ~LatticeRegular2d() {}

        /* virtual functions */
        virtual bool isSpanning(size_t = 1);

    protected:
        size_t b; // The boundary of any dimension of the lattice

};

#endif
