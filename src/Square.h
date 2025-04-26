/* Square.h
 * 08/01/05 Created this as a subclass of LatticeRegular2d.
*/

#ifndef SQUARE_H
#define SQUARE_H

#include "LatticeRegular2d.h"

class Square: public LatticeRegular2d
{
    public:
         Square(size_t len, size_t id = 0, short pbc = 2) : LatticeRegular2d( len, id ) 
         {
             Setup(pbc);
         }
         virtual ~Square() {}

    protected:
         short PBC;
         void Setup( short pbc = 2);
         void setNbrs0(size_t); /* no PBC */
         void setNbrsY(size_t); /* PBC in y-direction */
         void setNbrsXY(size_t); /* PBC in both directions */

};

#endif
