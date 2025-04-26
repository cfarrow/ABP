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

#include "LatticeRegular.h"

class LatticeRegular3d: public LatticeRegular
{
    public:
         LatticeRegular3d(size_t len, size_t id = 0) : LatticeRegular( len*len*len, id ) {}
         virtual ~LatticeRegular3d() {}

         /* virtual functions */
         virtual bool isSpanning(size_t = 1);

};
#endif
