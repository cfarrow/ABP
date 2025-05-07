/* LatticeRegular3d.h
 * This class is meant to be a superclass for 3D regular graphs.
 * It's neighbor retrieval scheme does not work with random graphs.
 * A proper subclass will only need to define the setNbrs() function,
 * which is inline. This means that there is no need for a cpp file.
 * The derived class may require additional items in the constructor.
 * See Cubic.h for an example.
 *
 * 08/01/05 Created this as a subclass of Lattice.
*/

#ifndef LATTICEREGULAR3D_H
#define LATTICEREGULAR3D_H

#include <cmath>
#include "LatticeRegular.h"


class LatticeRegular3d: public LatticeRegular
{
    public:
        LatticeRegular3d(size_t len, size_t id=0, size_t nn=0) 
        : LatticeRegular(len*len*len, id, nn) {
            length = static_cast<size_t>(ceil(pow(num_sites, 1.0/3.0)));
            b = length - 1;
        }
        virtual ~LatticeRegular3d() {}

        /* virtual functions */
        virtual bool isSpanning(size_t = 1);

    protected:
        size_t b; // The boundary of any dimension of the lattice

};

#endif
