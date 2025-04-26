/* BCC.h
 * 03/06/07 Created this as a subclass of LatticeRegular3d.
*/

#ifndef BCC_H
#define BCC_H

#include "LatticeRegular3d.h"

class BCC: public LatticeRegular3d
{
    public:
         BCC(size_t len, size_t id = 0, short pbc = 3)
             : LatticeRegular3d( len, id ) 
         {
             Setup(pbc);
         }
         virtual ~BCC() {}

    protected:
         short PBC;
         void Setup(short pbc = 3);
         void setNbrs0(size_t); /* no PBC */
         void setNbrsZ(size_t); /* PBC in one direction */
         void setNbrsYZ(size_t); /* PBC in two directions */
         void setNbrsXYZ(size_t); /* PBC in all directions */

};
#endif
