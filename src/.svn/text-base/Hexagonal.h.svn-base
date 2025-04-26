/* Hexagonal.h
 * 01/30/07 Created this as a subclass of LatticeRegular2d.
*/

#ifndef HEXAGONAL_H
#define HEXAGONAL_H

#include "LatticeRegular2d.h"

class Hexagonal: public LatticeRegular2d
{
    public:
         Hexagonal(size_t len, size_t id = 0, short pbc = 2) : LatticeRegular2d( len, id ) 
         {
             Setup(pbc);
         }
         virtual ~Hexagonal() {}

    protected:
         short PBC;
         void Setup( short pbc = 2 );
         void setNbrs0(size_t); /* no PBC */
         void setNbrsY(size_t); /* PBC in y-direction */
         void setNbrsXY(size_t); /* PBC in both directions */

};
#endif
