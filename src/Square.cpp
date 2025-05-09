#include "Square.h"
#include "Point.h"


/* PBC in both directions */

template<>
size_t Square<2>::getNumNeighbors(size_t) {
    return 4; 
}


template<>
size_t Square<2>::setNbrs(size_t i) {
    Point2d p{i, length};
    nbrs[1] = p.shift(-1,  0);
    nbrs[0] = p.shift( 1,  0);
    nbrs[3] = p.shift( 0, -1);
    nbrs[2] = p.shift( 0,  1);
    return 4;
}

/* PBC in one direction:
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1)   */

template<>
size_t Square<1>::getNumNeighbors(size_t i) {
    Point2d p{i, length};
    size_t nn = 4;
    if(p.x == 0 || p.x == b) nn -= 1;
    return nn;
}


template<>
size_t Square<1>::setNbrs(size_t i) {
    /* Get the coordinates of the lattice site */
    Point2d p{i, length};
    size_t n{0};

    if(p.x != 0) nbrs[n++] = p.shift(-1,  0);
    if(p.x != b) nbrs[n++] = p.shift( 1,  0);
    nbrs[n++] = p.shift( 0,  1);
    nbrs[n++] = p.shift( 0, -1);
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
size_t Square<0>::getNumNeighbors(size_t i) {
    Point2d p{i, length};
    size_t nn = 4;
    if(p.x == 0 || p.x == b) nn -= 1;
    if(p.y == 0 || p.y == b) nn -= 1;
    return nn;
}


template<>
size_t Square<0>::setNbrs(size_t i) {
    Point2d p{i, length};
    size_t n{0};

    if(p.x != 0) nbrs[n++] = p.shift(-1,  0);
    if(p.x != b) nbrs[n++] = p.shift( 1,  0);
    if(p.y != 0) nbrs[n++] = p.shift( 0, -1);
    if(p.y != b) nbrs[n++] = p.shift( 0, -1);
    return n;
}


// Let the compiler know we want these
template class Square<0>;
template class Square<1>;
template class Square<2>;