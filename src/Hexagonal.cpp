/* Hexagonal.h
 * 01/30/07 Created this as a subclass of LatticeRegular2d.
*/

#include "Hexagonal.h"
#include "Point.h"


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
    Point2d p{i, length};
    bool left = (p.x + p.y) % 2;
    nbrs[0] = p.shift(0, -1);
    nbrs[1] = p.shift(0,  1);
    nbrs[2] = p.shift(left ? -1 : 1, 0);
}


/* PBC in the Y direction */

template<>
size_t Hexagonal<1>::getNumNeighbors(size_t i) {
    Point2d p{i, length};
    bool left = (p.x + p.y) % 2;
    size_t nn = 3;
    if( left && p.x == 0)   nn -= 1;
    if(!left && p.x == b)   nn -= 1;
    return nn;
}


template<>
void Hexagonal<1>::setNbrs(size_t i) 
{
    Point2d p{i, length};
    bool left = (p.x + p.y) % 2;
    size_t n{0};

    nbrs[n++] = p.shift( 0, -1);
    nbrs[n++] = p.shift( 0,  1);
    if( left && p.x != 0)   nbrs[n++] = p.shift(-1, 0);
    if(!left && p.x != b)   nbrs[n++] = p.shift( 1, 0);
}


/* PBC in no directions: */

template<>
size_t Hexagonal<0>::getNumNeighbors(size_t i) {
    Point2d p{i, length};
    bool left = (p.x + p.y) % 2;
    size_t nn = 3;
    if( left && p.x == 0)       nn -= 1;
    if(!left && p.x == b)       nn -= 1;
    if(p.y == 0 || p.y == b)    nn -= 1;
    return nn;
}


template<>
void Hexagonal<0>::setNbrs(size_t i) 
{
    Point2d p{i, length};
    size_t n{0};
    bool left = (p.x + p.y) % 2;

    if(p.y != 0)            nbrs[n++] = p.shift( 0, -1);
    if(p.y != b)            nbrs[n++] = p.shift( 0,  1);
    if( left && p.x != 0)   nbrs[n++] = p.shift(-1,  0);
    if(!left && p.x != b)   nbrs[n++] = p.shift( 1,  0);
}


// Let the compiler know we want these
template class Hexagonal<0>;
template class Hexagonal<1>;
template class Hexagonal<2>;
