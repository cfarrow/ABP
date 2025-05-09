#include "UJack.h"
#include "Point.h"


template<>
size_t UJack<2>::getNumNeighbors(size_t) {
    return 8; 
}


template<>
size_t UJack<2>::setNbrs(size_t i) 
{
    Point2d p{i, length};
    nbrs[0] = p.shift(-1,  0);
    nbrs[1] = p.shift( 1,  0);
    nbrs[2] = p.shift( 0, -1);
    nbrs[3] = p.shift( 0,  1);
    nbrs[4] = p.shift(-1,  1);
    nbrs[5] = p.shift( 1, -1);
    nbrs[6] = p.shift( 1,  1);
    nbrs[7] = p.shift(-1, -1);
    return 8;
}


/* PBC in one direction:
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1)   */
template<>
size_t UJack<1>::getNumNeighbors(size_t i) {
    Point2d p{i, length};
    size_t nn = 8;
    if(p.x == 0 || p.x == b) nn -= 3;
    return nn;
}


template<>
size_t UJack<1>::setNbrs(size_t i) 
{
    Point2d p{i, length};
    size_t n{0};

    if(p.x != 0){
        nbrs[n++] = p.shift(-1,  0);
        nbrs[n++] = p.shift(-1,  1);
        nbrs[n++] = p.shift(-1, -1);
    }
    if(p.x != b){
        nbrs[n++] = p.shift( 1,  0);
        nbrs[n++] = p.shift( 1, -1);
        nbrs[n++] = p.shift( 1,  1);
    }
    nbrs[n++] = p.shift( 0, -1);
    nbrs[n++] = p.shift( 0,  1);
    return n;
}


/* PBC in no directions:
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1)   
 * We must also check whether if the y-value is 0 or length-1, etc.
 * It is helpful to visualize the lattice extending from [0,0] to
 * the right and up.
 * */
template<>
size_t UJack<0>::getNumNeighbors(size_t i) {
    Point2d p{i, length};
    size_t nn = 8;
    bool x_boundary = (p.x == 0 || p.x == b);
    bool y_boundary = (p.y == 0 || p.y == b);
    if(x_boundary) nn -= 3;
    if(y_boundary) nn -= 3;
    // If we're in a corner, then we counted a diagonal neighbor twice
    if(x_boundary && y_boundary) nn += 1;
    return nn;
}


template<>
size_t UJack<0>::setNbrs(size_t i) 
{
    Point2d p{i, length};
    size_t n{0};
    if(p.x != 0)                nbrs[n++] = p.shift(-1,  0);
    if(p.x != b)                nbrs[n++] = p.shift( 1,  0);
    if(p.y != 0)                nbrs[n++] = p.shift( 0, -1);
    if(p.y != b)                nbrs[n++] = p.shift( 0,  1);
    if(p.x != 0 && p.y != b)    nbrs[n++] = p.shift(-1,  1);
    if(p.x != b && p.x != 0)    nbrs[n++] = p.shift( 1, -1);
    if(p.x != b && p.y != b)    nbrs[n++] = p.shift( 1,  1);
    if(p.x != 0 && p.y != 0)    nbrs[n++] = p.shift(-1, -1);
    return n;
}


// Let the compiler know we want these
template class UJack<0>;
template class UJack<1>;
template class UJack<2>;