#include "Cubic.h"
#include "Point.h"

/* PBC in all directions */

template<>
size_t Cubic<3>::getNumNeighbors(size_t) {
    return 6;
}


template <>
size_t Cubic<3>::setNbrs(size_t i) 
{
    Point3d p{i, length};
    nbrs[0] = p.shift(-1,  0,  0);
    nbrs[1] = p.shift( 1,  0,  0);
    nbrs[2] = p.shift( 0, -1,  0);
    nbrs[3] = p.shift( 0,  1,  0);
    nbrs[4] = p.shift( 0,  0, -1);
    nbrs[5] = p.shift( 0,  0,  1);
    return 6;
}


/* PBC in two directions: 
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1).   
 * This results in PBC in the Y and Z directions.
 */
template<>
size_t Cubic<2>::getNumNeighbors(size_t i) {
    Point3d p{i, length};
    size_t nn = 6;
    if(p.x == 0 || p.x == length-1) nn -= 1;
    return nn;
}


template <>
size_t Cubic<2>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    Point3d p{i, length};
    size_t n{0};
    if(p.x != 0)            nbrs[n++] = p.shift(-1,  0,  0);
    if(p.x != length-1)     nbrs[n++] = p.shift( 1,  0,  0);
    nbrs[n++] = p.shift( 0, -1,  0);
    nbrs[n++] = p.shift( 0,  1,  0);
    nbrs[n++] = p.shift( 0,  0, -1);
    nbrs[n++] = p.shift( 0,  0,  1);
    return n;
}


/* PBC in one direction: 
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1).   
 * We do the same for y.
 * This results in PBC in the Z direction.
 */
template<>
size_t Cubic<1>::getNumNeighbors(size_t i) {
    Point3d p{i, length};
    size_t b = length - 1;
    size_t nn = 6;
    if(p.x == 0 || p.x == b) nn -= 1;
    if(p.y == 0 || p.y == b) nn -= 1;
    return nn;
}


template <>
size_t Cubic<1>::setNbrs(size_t i) 
{
    Point3d p{i, length};
    size_t n{0};
    if(p.x != 0) nbrs[n++] = p.shift(-1,  0,  0);
    if(p.x != b) nbrs[n++] = p.shift( 1,  0,  0);
    if(p.y != 0) nbrs[n++] = p.shift( 0, -1,  0);
    if(p.y != b) nbrs[n++] = p.shift( 0,  1,  0);
    nbrs[n++] = p.shift( 0,  0, -1);
    nbrs[n++] = p.shift( 0,  0,  1);
    return n;
}


/* PBC in no direction: 
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1).   
 * We do the same for y and z.
 */
template<>
size_t Cubic<0>::getNumNeighbors(size_t i) {
    Point3d p{i, length};
    size_t b = length - 1;
    size_t nn = 6;
    if(p.x == 0 || p.x == b) nn -= 1;
    if(p.y == 0 || p.y == b) nn -= 1;
    if(p.z == 0 || p.z == b) nn -= 1;
    return nn;
}


template <>
size_t Cubic<0>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    Point3d p{i, length};
    size_t n{0};
    if(p.x != 0) nbrs[n++] = p.shift(-1,  0,  0);
    if(p.x != b) nbrs[n++] = p.shift( 1,  0,  0);
    if(p.y != 0) nbrs[n++] = p.shift( 0, -1,  0);
    if(p.y != b) nbrs[n++] = p.shift( 0,  1,  0);
    if(p.z != 0) nbrs[n++] = p.shift( 0,  0, -1);
    if(p.z != b) nbrs[n++] = p.shift( 0,  0,  1);
    return n;
}


// Let the compiler know we want these
template class Cubic<0>;
template class Cubic<1>;
template class Cubic<2>;
template class Cubic<3>;