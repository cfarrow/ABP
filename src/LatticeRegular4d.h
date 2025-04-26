/* LatticeRegular4d.h
 *
 * 03/15/07 Created this as a subclass of LatticeRegular.
*/

#ifndef LATTICEREGULAR4D_H
#define LATTICEREGULAR4D_H

#include "LatticeRegular.h"

class LatticeRegular4d: public LatticeRegular
{
    public:
         LatticeRegular4d(size_t len, size_t id = 0) : LatticeRegular( len*len*len*len, id ) {}
         virtual ~LatticeRegular4d() {}

         /* virtual functions */
         virtual bool isSpanning(size_t = 1);

};
#endif
