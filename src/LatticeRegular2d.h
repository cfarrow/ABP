/* LatticeRegular2d.h
 * This class is meant to be a superclass for 2D regular graphs.
 * It's neighbor retrieval scheme does not work with random graphs.
 * A proper subclass will only need to define the setNbrs() function,
 * which is inline. This means that there is no need for a cpp file.
 * The derived class may require additional items in the constructor.
 * See Triangular.h for an example.
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
        LatticeRegular2d(size_t len, size_t id = 0, size_t nn = 0)
        : LatticeRegular(len*len, id) {
           length = static_cast< size_t >(ceil(sqrt(num_sites)));
           num_neighbors = nn;
           nbrs.resize(num_neighbors);
           last = num_sites;

        }
        virtual ~LatticeRegular2d() {}

        /* virtual functions */
        virtual bool isSpanning(size_t = 1);

};

#endif
