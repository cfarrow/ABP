#include "Square.h"
#include <vector>
#include <cmath>

void Square::Setup(short pbc)
{
    /* If num_sites was changed in the constructor, this makes sure that
     * length is properly defined.
     */
    length = static_cast< size_t >( ceil(sqrt(num_sites)) );
    num_neighbors = 4;
    nbrs.resize(num_neighbors);
    last = num_sites;

    // Set the function that will return the neighbors
    PBC = pbc;
    if( PBC >= 2 ) setNbrs = (LRFptr) &Square::setNbrsXY;
    else if( PBC == 1) setNbrs =  (LRFptr) &Square::setNbrsY;
    else setNbrs = (LRFptr) &Square::setNbrs0;
}


/* PBC in both directions */
void Square::setNbrsXY(size_t i) 
{
	/* populates the neighbor array */
	size_t x_val = 0, y_val = 0;

    /* Get the coordinates of the lattice site */
    x_val = i % length;
    y_val = ( i - x_val ) / length;

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
void Square::setNbrsY(size_t i) 
{
	/* populates the neighbor array */
	size_t x_val = 0, y_val = 0;

    /* Get the coordinates of the lattice site */
    x_val = i % length;
    y_val = ( i - x_val ) / length;

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
void Square::setNbrs0(size_t i) 
{
	/* populates the neighbor array */
	size_t x_val = 0, y_val = 0;

    /* Get the coordinates of the lattice site */
    x_val = i % length;
    y_val = ( i - x_val ) / length;

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
