#include "UJack.h"


/* For regular graphs, the convention is to create an array of neighbors
 * of site i with this function. It then is called whenever a neighbor
 * is needed. This saves space since a bond array is not required. This
 * will not work with random graphs, however.
 * THIS IS NOT DEFINED IN LATTICE_H SINCE IT IS SPECIFIC TO REGULAR GRAPHS.
 */
template<>
size_t UJack<2>::getNumNeighbors(size_t) {
    return 8; 
}


template<>
void UJack<2>::setNbrs(size_t i) 
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
		
    /* up and right, x+1 & y+1 */
    nbrs[6] = ( x_val + 1 )%length + (( y_val + 1 )%length) * length;

    /* down and left, x-1 & y-1 */
    nbrs[7] = ( x_val - 1 - length )%length + (( y_val - 1 + length )%length) * length;
}


/* PBC in one direction:
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1)   */
template<>
size_t UJack<1>::getNumNeighbors(size_t i) {
	size_t x_val = i % length;
    size_t nn = 8;
    if(x_val == 0 || x_val == length-1) nn -= 3;
    return nn;
}


template<>
void UJack<1>::setNbrs(size_t i) 
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

    /* up and right, x+1 & y+1 */
    if( x_val != length -1 ) {
        nbrs.push_back( ( x_val + 1 ) + (( y_val + 1 )%length) * length );
    }

    /* down and left, x-1 & y-1 */
    if( x_val != 0 ) {
        nbrs.push_back( ( x_val - 1 ) + (( y_val - 1 + length )%length) * length );
    }
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
	size_t x_val = i % length;
    size_t y_val = ( i - x_val ) / length;
    size_t nn = 8;

    bool x_boundary = (x_val == 0 || x_val == length-1);
    bool y_boundary = (y_val == 0 || y_val == length-1);

    if(x_boundary) nn -= 3;
    if(y_boundary) nn -= 3;
    // If we're in a corner, then we counted a diagonal neighbor twice
    if(x_boundary && y_boundary) nn += 1;
    return nn;
}


template<>
void UJack<0>::setNbrs(size_t i) 
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

    /* up and right, x+1 & y+1 */
    if( x_val != length -1  && y_val != length -1) {
        nbrs.push_back( ( x_val + 1 ) + ( y_val + 1 ) * length );
    }

    /* down and left, x-1 & y-1 */
    if( x_val != 0  && y_val != 0 ) {
        nbrs.push_back( ( x_val - 1 ) + ( y_val - 1 ) * length );
    }
}


// Let the compiler know we want these
template class UJack<0>;
template class UJack<1>;
template class UJack<2>;