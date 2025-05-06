#include "Triangular.h"
#include "Point.h"


/* PBC in both directions. Every site is treated the same */

template<>
size_t Triangular<2>::getNumNeighbors(size_t) {
    return 6; 
}


template<>
void Triangular<2>::setNbrs(size_t i) 
{
    Point2d p{i, length};
    nbrs[0] = p.shift( 1,  0);
    nbrs[1] = p.shift( 0,  1);
    nbrs[2] = p.shift(-1,  0);
    nbrs[3] = p.shift( 0, -1);
    nbrs[4] = p.shift(-1,  1);
    nbrs[5] = p.shift( 1, -1);
}


/* PBC in one direction:
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1)   
 * This results in PBC in the Y direction.
 */

template<>
size_t Triangular<1>::getNumNeighbors(size_t i) {
    Point2d p{i, length};
    size_t nn = 6;
    if(p.x == 0 || p.x == b) nn -= 2;
    return nn;
}


template<>
void Triangular<1>::setNbrs(size_t i) 
{
    Point2d p{i, length};
    nbrs.clear();
    nbrs.push_back(p.shift( 0,  1));
    nbrs.push_back(p.shift( 0, -1));
    if(p.x != 0) {
        nbrs.push_back(p.shift(-1,  0));
        nbrs.push_back(p.shift(-1,  1));
    }
    if(p.x != b) {
        nbrs.push_back(p.shift( 1,  0));
        nbrs.push_back(p.shift( 1, -1));
    }
}


/* PBC in no directions:
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1)   
 * We must also check whether if the y-value is 0 or length-1, etc.
 */

template<>
size_t Triangular<0>::getNumNeighbors(size_t i) {
    Point2d p{i, length};
    size_t nn = 6;
    if(p.x == 0 || p.x == b) nn -= 2;
    if(p.y == 0 || p.y == b) nn -= 2;
    
	// The upper-left and lower-right corners are convex, so +1 nbr
    if(p.x == 0 && p.y == b) nn += 1;
    if(p.x == b && p.y == 0) nn += 1;
    return nn;
}


template<>
void Triangular<0>::setNbrs(size_t i) 
{
    Point2d p{i, length};
    bool left = p.x == 0;
    bool right = p.x == b;
    bool botm = p.y == 0;
    bool top = p.y == b;

    nbrs.clear();
    if(!left)               nbrs.push_back(p.shift(-1,  0));
    if(!right)              nbrs.push_back(p.shift( 1,  0));
    if(!botm)               nbrs.push_back(p.shift( 0, -1));
    if(!top)                nbrs.push_back(p.shift( 0,  1));
    if(!left && !top)       nbrs.push_back(p.shift(-1,  1));
    if(!right && !botm)     nbrs.push_back(p.shift( 1, -1));
}


// Let the compiler know we want these
template class Triangular<0>;
template class Triangular<1>;
template class Triangular<2>;