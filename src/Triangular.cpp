#include "Triangular.h"


/* PBC in both directions. Every site is treated the same */

template<>
size_t Triangular<2>::getNumNeighbors(size_t) {
    return 6; 
}


template<>
void Triangular<2>::setNbrs(size_t i) 
{

    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( i - x_val ) / length;

    /* right, x + 1 */
    nbrs[0] = ( x_val + 1 ) % length + y_val * length;
    
    /* up, y + 1 */
    nbrs[1] = x_val + (( y_val + 1 )%length) * length;

    /* left, x - 1 */
    nbrs[2] = ( x_val - 1 + length )%length + y_val * length;
    
    /* down, y - 1 */
    nbrs[3] = x_val + (( y_val - 1 + length )%length) * length;

	/* up and left, x-1 & y+1 */
	nbrs[4] = ( x_val - 1 + length )%length + (( y_val + 1 )%length) * length;

	/* down and right, x+1 & y-1 */
	nbrs[5] = ( x_val + 1 )%length + (( y_val - 1 + length )%length) * length;
}


/* PBC in one direction:
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1)   
 * This results in PBC in the Y direction.
 */

template<>
size_t Triangular<1>::getNumNeighbors(size_t i) {
	size_t x_val = i % length;
    size_t nn = 6;
    if(x_val == 0 || x_val == length-1) nn -= 2;
    return nn;
}


template<>
void Triangular<1>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( i - x_val ) / length;

    nbrs.clear();

    /* right, x + 1 */
    if( x_val != length -1 ) {
        nbrs.push_back( ( x_val + 1 ) + y_val * length );
    }
    
    /* up, y + 1 */
    nbrs.push_back( x_val + (( y_val + 1 )%length) * length );

    /* left, x - 1 */
    if( x_val != 0 ) {
        nbrs.push_back( ( x_val - 1 ) + y_val * length);
    }
    
    /* down, y - 1 */
    nbrs.push_back( x_val + (( y_val - 1 + length )%length) * length);

	/* up and left, x-1 & y+1 */
    if( x_val != 0 ) {
        nbrs.push_back( ( x_val - 1 ) + (( y_val + 1 )%length) * length);
    }

	/* down and right, x+1 & y-1 */
    if( x_val != length -1 ) {
        nbrs.push_back( ( x_val + 1 ) + (( y_val - 1 + length )%length) * length);
    }
}


/* PBC in no directions:
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1)   
 * We must also check whether if the y-value is 0 or length-1, etc.
 */

template<>
size_t Triangular<0>::getNumNeighbors(size_t i) {
    size_t x_val = i % length;
    size_t y_val = ( i - x_val ) / length;

    size_t nn = 6;
    if(x_val == 0 || x_val == length-1) nn -= 2;
    if(y_val == 0 || y_val == length-1) nn -= 2;
    
	// The upper-left and lower-right corners are convex, so +1 nbr
    if( x_val == 0 && y_val == length-1) {
        nn += 1;
    }
    if( x_val == length -1 && y_val == 0) {
        nn += 1;
    }
    return nn;
}

template<>
void Triangular<0>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( i - x_val ) / length;

    nbrs.clear();

    /* right, x + 1 */
    if( x_val != length -1 ) {
        nbrs.push_back( ( x_val + 1 ) + y_val * length );
    }
    
    /* up, y + 1 */
    if( y_val != length -1 ) {
        nbrs.push_back( x_val + ( y_val + 1 ) * length );
    }

    /* left, x - 1 */
    if( x_val != 0 ) {
        nbrs.push_back( ( x_val - 1 ) + y_val * length);
    }
    
    /* down, y - 1 */
    if( y_val != 0 ) {
        nbrs.push_back( x_val + ( y_val - 1 ) * length);
    }

	/* up and left, x-1 & y+1 */
    if( x_val != 0 && y_val != length-1) {
        nbrs.push_back( ( x_val - 1 ) + (( y_val + 1 )%length) * length);
    }

	/* down and right, x+1 & y-1 */
    if( x_val != length -1 && y_val != 0) {
        nbrs.push_back( ( x_val + 1 ) + (( y_val - 1 + length )%length) * length);
    }
}


// Let the compiler know we want these
template class Triangular<0>;
template class Triangular<1>;
template class Triangular<2>;