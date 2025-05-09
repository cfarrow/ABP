#include "Cubic4d.h"
#include "Point.h"

/* PBC in all directions */

template<>
size_t Cubic4d<4>::getNumNeighbors(size_t) {
    return 8;
}


template<>
size_t Cubic4d<4>::setNbrs(size_t i) 
{
    Point4d p{i, length};
    nbrs[0] = p.shift(-1,  0,  0,  0);
    nbrs[1] = p.shift( 1,  0,  0,  0);
    nbrs[2] = p.shift( 0, -1,  0,  0);
    nbrs[3] = p.shift( 0,  1,  0,  0);
    nbrs[4] = p.shift( 0,  0, -1,  0);
    nbrs[5] = p.shift( 0,  0,  1,  0);
    nbrs[6] = p.shift( 0,  0,  0, -1);
    nbrs[7] = p.shift( 0,  0,  0,  1);
    return 8;
}


/* PBC in three directions: 
 * We check whether the w-value is 0 or length-1. If it is then we don't connect
 * its neighbor across the boundary.
 */
template<>
size_t Cubic4d<3>::getNumNeighbors(size_t i) {
    Point4d p{i, length};
    size_t nn = 8;
    if(p.w == 0 || p.w == b) nn -= 1;
    return nn;
}


template<>
size_t Cubic4d<3>::setNbrs(size_t i) 
{
    Point4d p{i, length};
    size_t n{0};
    if(p.w != 0)    nbrs[n++] = p.shift(-1,  0,  0,  0);
    if(p.w != b)    nbrs[n++] = p.shift( 1,  0,  0,  0);
    nbrs[n++] = p.shift( 0, -1,  0,  0);
    nbrs[n++] = p.shift( 0,  1,  0,  0);
    nbrs[n++] = p.shift( 0,  0, -1,  0);
    nbrs[n++] = p.shift( 0,  0,  1,  0);
    nbrs[n++] = p.shift( 0,  0,  0, -1);
    nbrs[n++] = p.shift( 0,  0,  0,  1);
    return n;
}


/* PBC in two directions: 
 * We check whether the w- or x-value is 0 or length-1. If it is then we don't
 * connect its neighbor across the boundary.
 */
template<>
size_t Cubic4d<2>::getNumNeighbors(size_t i) {
    Point4d p{i, length};
    size_t nn = 8;
    if(p.w == 0 || p.w == b) nn -= 1;
    if(p.x == 0 || p.x == b) nn -= 1;
    return nn;
}


template<>
size_t Cubic4d<2>::setNbrs(size_t i) 
{
    Point4d p{i, length};
    size_t n{0};
    if(p.w != 0)    nbrs[n++] = p.shift(-1,  0,  0,  0);
    if(p.w != b)    nbrs[n++] = p.shift( 1,  0,  0,  0);
    if(p.x != 0)    nbrs[n++] = p.shift( 0, -1,  0,  0);
    if(p.x != b)    nbrs[n++] = p.shift( 0,  1,  0,  0);
    nbrs[n++] = p.shift( 0,  0, -1,  0);
    nbrs[n++] = p.shift( 0,  0,  1,  0);
    nbrs[n++] = p.shift( 0,  0,  0, -1);
    nbrs[n++] = p.shift( 0,  0,  0,  1);
    return n;
}

/* PBC in one directions: 
 * We check whether the w-, x-, or y-value is 0 or length-1. If it is then we
 * don't connect its neighbor across the boundary.
 */
template<>
size_t Cubic4d<1>::getNumNeighbors(size_t i) {
    Point4d p{i, length};
    size_t nn = 8;
    if(p.w == 0 || p.w == b) nn -= 1;
    if(p.x == 0 || p.x == b) nn -= 1;
    if(p.y == 0 || p.y == b) nn -= 1;
    return nn;
}

template<>
size_t Cubic4d<1>::setNbrs(size_t i) 
{
    Point4d p{i, length};
    size_t n{0};
    if(p.w != 0)    nbrs[n++] = p.shift(-1,  0,  0,  0);
    if(p.w != b)    nbrs[n++] = p.shift( 1,  0,  0,  0);
    if(p.x != 0)    nbrs[n++] = p.shift( 0, -1,  0,  0);
    if(p.x != b)    nbrs[n++] = p.shift( 0,  1,  0,  0);
    if(p.y != 0)    nbrs[n++] = p.shift( 0,  0, -1,  0);
    if(p.y != b)    nbrs[n++] = p.shift( 0,  0,  1,  0);
    nbrs[n++] = p.shift( 0,  0,  0, -1);
    nbrs[n++] = p.shift( 0,  0,  0,  1);
    return n;
}


/* PBC in no directions: 
 * We check whether the w-, x-, y-, or z-value is 0 or length-1. If it is then
 * we don't connect its neighbor across the boundary.
 */
template<>
size_t Cubic4d<0>::getNumNeighbors(size_t i) {
    Point4d p{i, length};
    size_t nn = 8;
    if(p.w == 0 || p.w == b) nn -= 1;
    if(p.x == 0 || p.x == b) nn -= 1;
    if(p.y == 0 || p.y == b) nn -= 1;
    if(p.z == 0 || p.z == b) nn -= 1;
    return nn;
}


template<>
size_t Cubic4d<0>::setNbrs(size_t i) 
{
    Point4d p{i, length};
    size_t n{0};
    if(p.w != 0)    nbrs[n++] = p.shift(-1,  0,  0,  0);
    if(p.w != b)    nbrs[n++] = p.shift( 1,  0,  0,  0);
    if(p.x != 0)    nbrs[n++] = p.shift( 0, -1,  0,  0);
    if(p.x != b)    nbrs[n++] = p.shift( 0,  1,  0,  0);
    if(p.y != 0)    nbrs[n++] = p.shift( 0,  0, -1,  0);
    if(p.y != b)    nbrs[n++] = p.shift( 0,  0,  1,  0);
    if(p.z != 0)    nbrs[n++] = p.shift( 0,  0,  0, -1);
    if(p.z != b)    nbrs[n++] = p.shift( 0,  0,  0,  1);
    return n;
}


// Let the compiler know we want these
template class Cubic4d<0>;
template class Cubic4d<1>;
template class Cubic4d<2>;
template class Cubic4d<3>;