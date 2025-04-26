/* Cubic.h
 * 08/01/05 Created this as a subclass of LatticeRegular3d.
*/

#ifndef CUBIC_H
#define CUBIC_H

#include "LatticeRegular3d.h"

class Cubic: public LatticeRegular3d
{
    public:
         Cubic(size_t len, size_t id = 0, short pbc = 3) 
             : LatticeRegular3d( len, id ) 
         {
             Setup(pbc);
         }
         virtual ~Cubic() {}

    protected:
         short PBC;
         void Setup(short pbc = 3);
         void setNbrs0(size_t); /* no PBC */
         void setNbrsZ(size_t); /* PBC in one direction */
         void setNbrsYZ(size_t); /* PBC in two directions */
         void setNbrsXYZ(size_t); /* PBC in all directions */

};
#endif
