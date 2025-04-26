/* Hexagonal.h
 * 01/30/07 Created this as a subclass of LatticeRegular2d.
*/

#include "Hexagonal.h"
#include <vector>
#include <cmath>

void Hexagonal::Setup(short pbc)
{
    /* If num_sites was changed in the constructor, this makes sure that
     * length is properly defined.
     */
    length = static_cast< size_t >( ceil(sqrt(num_sites)) );
    num_neighbors = 3;
    nbrs.resize(num_neighbors);
    last = num_sites;

    PBC = pbc;
    if( PBC >= 2 ) setNbrs = (LRFptr) &Hexagonal::setNbrsXY;
    else if( PBC ==1 )setNbrs =  (LRFptr) &Hexagonal::setNbrsY;
    else setNbrs = (LRFptr) &Hexagonal::setNbrs0;
}

/* PBC in both directions */
void Hexagonal::setNbrsXY(size_t i) 
{
	/* populates the neighbor array */
	size_t x_val = 0, y_val = 0;

    /* Get the coordinates of the lattice site */
    x_val = i % length;
    y_val = ( i - x_val ) / length;

    /* up, y + 1 */
    nbrs[0] = x_val + (( y_val + 1 )%length) * length;

    /* down, y - 1 */
    nbrs[1] = x_val + (( y_val - 1 + length )%length) * length;

    /* left/right, x+/- 1
     * if x_val + y_val is even, the site is connected to the right
     * if it is odd, the site is connected to the left
     */
    nbrs[2] = (length + x_val + ((x_val+y_val)%2 ? 1 : -1))%length + y_val * length;
}

/* PBC in the Y-direction direction: */
void Hexagonal::setNbrsY(size_t i) 
{
	/* populates the neighbor array */
	size_t x_val = 0, y_val = 0;

    /* Get the coordinates of the lattice site */
    x_val = i % length;
    y_val = ( i - x_val ) / length;

    nbrs.clear();

    /* up, y + 1 */
    nbrs.push_back( x_val + (( y_val + 1 )%length) * length );

    /* down, y - 1 */
    nbrs.push_back( x_val + (( y_val - 1 + length )%length) * length );

    /* left, x - 1 */
    if(!((x_val+y_val)%2) and x_val != 0) {
    nbrs.push_back(  x_val - 1 + y_val * length );
    }

    /* right, x + 1 */
    if((x_val+y_val)%2 and x_val != length-1) {
    nbrs.push_back( x_val + 1 + y_val * length );
    }
}

/* PBC in no directions: */
void Hexagonal::setNbrs0(size_t i) 
{
	/* populates the neighbor array */
	size_t x_val = 0, y_val = 0;

    /* Get the coordinates of the lattice site */
    x_val = i % length;
    y_val = ( i - x_val ) / length;

    nbrs.clear();

    /* up, y + 1 */
    if(y_val != length-1) {
        nbrs.push_back( x_val + ( y_val + 1 ) * length );
    }

    /* down, y - 1 */
    if(y_val != 0) {
        nbrs.push_back( x_val + ( y_val - 1) * length );
    }

    /* left, x - 1 */
    if(!((x_val+y_val)%2) and x_val != 0) {
    nbrs.push_back( x_val - 1 + y_val * length );
    }

    /* left, x + 1 */
    if((x_val+y_val)%2 and x_val != length-1) {
    nbrs.push_back( x_val + 1 + y_val * length );
    }
}
