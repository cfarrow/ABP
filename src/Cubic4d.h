/* Cubic4d.h
 * 03/15/07 Created this as a subclass of LatticeRegular4d.
*/

#ifndef CUBIC4D_H
#define CUBIC4D_H

#include "LatticeRegular4d.h"

class Cubic4d: public LatticeRegular4d
{
    public:
         Cubic4d(size_t len, size_t id = 0, short pbc = 4) 
             : LatticeRegular4d( len, id ) 
         {
             Setup(pbc);
         }
         virtual ~Cubic4d() {}

    protected:
         short PBC;
         void Setup(short pbc = 4);
         void setNbrs0(size_t); /* no PBC */
         void setNbrsZ(size_t); /* PBC in one direction */
         void setNbrsYZ(size_t); /* PBC in two directions */
         void setNbrsXYZ(size_t); /* PBC in three directions */
         void setNbrsWXYZ(size_t); /* PBC in all directions */

};
#endif
