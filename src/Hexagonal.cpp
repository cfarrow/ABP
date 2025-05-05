/* Hexagonal.h
 * 01/30/07 Created this as a subclass of LatticeRegular2d.
*/

#include "Hexagonal.h"
#include <vector>
#include <cmath>


/* PBC in both direction

A hexagonal lattice has three neighbors: 
- one directly above (y + 1)
- one directly below (y - 1)
- one to the right or left (right if x + y is even, left otherwise)

*/

template<>
size_t Hexagonal<2>::getNumNeighbors(size_t) {
    return 3; 
}


template<>
void Hexagonal<2>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( i - x_val ) / length;

    /* up, y + 1 */
    nbrs[0] = x_val + (( y_val + 1 )%length) * length;

    /* down, y - 1 */
    nbrs[1] = x_val + (( y_val - 1 + length )%length) * length;

    /* left/right, x+/- 1 */
    nbrs[2] = (length + x_val + ((x_val+y_val)%2 ? 1 : -1))%length + y_val * length;
}


/* PBC in the Y-direction direction: */

template<>
size_t Hexagonal<1>::getNumNeighbors(size_t i) {
	size_t x_val = i % length;
    size_t y_val = ( i - x_val ) / length;
    size_t nn = 3;

    /* The upper and lower neighbors remain.
     * if x_val + y_val is even, the site is connected to the right
     * if it is odd, the site is connected to the left
     * the left boundary is even and the right boundary is odd
     * Thus, points on those boundaries are missing a neighbor if they are in an
     * odd row (y value).
    */
    if(y_val %2 == 1 && (x_val == 0 || x_val == length - 1)) nn -= 1;
    return nn;
}


template<>
void Hexagonal<1>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( i - x_val ) / length;

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

template<>
size_t Hexagonal<0>::getNumNeighbors(size_t i) {
	size_t x_val = i % length;
    size_t y_val = ( i - x_val ) / length;
    size_t nn = 3;

    /* if x_val + y_val is even, the site is connected to the right
     * if it is odd, the site is connected to the left
     * the left boundary is even and the right boundary is odd
     * Thus, points on those boundaries are missing a neighbor if they are in an
     * odd row (y value).
     * Points on the top and bottom boundaries are missing a neighbor.
    */
    if(y_val %2 == 1 && (x_val == 0 || x_val == length - 1)) nn -= 1;
    if(y_val == 0 || y_val == length - 1) nn -= 1;
    return nn;
}


template<>
void Hexagonal<0>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( i - x_val ) / length;

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


// Let the compiler know we want these
template class Hexagonal<0>;
template class Hexagonal<1>;
template class Hexagonal<2>;
