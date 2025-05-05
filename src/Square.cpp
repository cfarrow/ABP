#include "Square.h"


/* PBC in both directions */

template<>
size_t Square<2>::getNumNeighbors(size_t) {
    return 4; 
}


template<>
void Square<2>::setNbrs(size_t i) {
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
}

/* PBC in one direction:
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1)   */

template<>
size_t Square<1>::getNumNeighbors(size_t i) {
	size_t x_val = i % length;
    size_t nn = 4;
    if(x_val == 0 || x_val == length-1) nn -= 1;
    return nn;
}


template<>
void Square<1>::setNbrs(size_t i) {
    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( i - x_val ) / length;

    nbrs.clear();

    /* right, x + 1 */
    if( x_val != length-1 ) {
        nbrs.push_back( ( x_val + 1 ) + y_val * length );
    }
    
    /* up, y + 1 */
    nbrs.push_back( x_val + (( y_val + 1 )%length) * length );

    /* left, x - 1 */
    if( x_val != 0 ) {
        nbrs.push_back( ( x_val - 1 ) + y_val * length );
    }
    
    /* down, y - 1 */
    nbrs.push_back( x_val + (( y_val - 1 + length )%length) * length );
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
    size_t x_val = i % length;
    size_t y_val = ( i - x_val ) / length;

    size_t nn = 4;
    if(x_val == 0 || x_val == length-1) nn -= 1;
    if(y_val == 0 || y_val == length-1) nn -= 1;
    return nn;
}


template<>
void Square<0>::setNbrs(size_t i) {
    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( i - x_val ) / length;

    nbrs.clear();

    /* right, x + 1 */
    if( x_val != length-1 ) {
        nbrs.push_back( ( x_val + 1 ) + y_val * length );
    }
    
    /* up, y + 1 */
    if( y_val != length-1 ) {
        nbrs.push_back( x_val + ( y_val + 1 ) * length );
    }

    /* left, x - 1 */
    if( x_val != 0 ) {
        nbrs.push_back( ( x_val - 1 ) + y_val * length );
    }
    
    /* down, y - 1 */
    if( y_val != 0 ) {
        nbrs.push_back( x_val + ( y_val - 1 ) * length );
    }
}


// Let the compiler know we want these
template class Square<0>;
template class Square<1>;
template class Square<2>;