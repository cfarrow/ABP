/* Triangular.h
 * 08/01/05 Created this as a subclass of LatticeRegular2d.
*/

#ifndef TRIANGULAR_H
#define TRIANGULAR_H

#include "LatticeRegular2d.h"

class Triangular: public LatticeRegular2d
{
    public:
         Triangular(size_t len, size_t id = 0, short pbc = 2) : LatticeRegular2d( len, id ) 
         {
             Setup(pbc);
         }
         virtual ~Triangular() {}

    protected:
         short PBC;
         void Setup(short pbc = 2);
         void setNbrs0(size_t); /* no PBC */
         void setNbrsY(size_t); /* PBC in y-direction */
         void setNbrsXY(size_t); /* PBC in both directions */

};

#endif
