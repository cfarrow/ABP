#include "Cubic.h"
#include <vector>
#include <cmath>

void Cubic::Setup( short pbc )
{
    /* If num_sites was changed in the constructor, this makes sure that
     * length is properly defined.
     */
    length = static_cast< size_t >( ceil( pow(num_sites, 1.0/3.0)) );
    num_neighbors = 6;
    nbrs.resize(num_neighbors);
    last = num_sites;

    /* Set the boundary conditions and set the correct setNbrs function */
    PBC = pbc;
    if( PBC >= 3 ) setNbrs = (LRFptr) &Cubic::setNbrsXYZ;
    else if( PBC == 2 ) setNbrs = (LRFptr) &Cubic::setNbrsYZ;
    else if( PBC == 1 ) setNbrs = (LRFptr) &Cubic::setNbrsZ;
    else setNbrs = (LRFptr) &Cubic::setNbrs0;
}


/* PBC in all directions */
void Cubic::setNbrsXYZ(size_t i) 
{
    /* If periodic boundary conditions are to be added, this is the place! */

	/* populates the neighbor array */
	size_t x_val, y_val, z_val;
	x_val = y_val = z_val = 0;

    /* Get the coordinates of the lattice site */
    x_val = i % length;
    y_val = ( (i - x_val)/length ) % length;
    z_val = (i - x_val - y_val*length) / ( length * length );

    /* right, x + 1 */
    nbrs[0] = ( x_val + 1 )%length + y_val * length + z_val * length * length;
    
    /* up, y + 1 */
    nbrs[1] = x_val + (( y_val + 1 )%length) * length + z_val * length * length;

    /* left, x - 1 */
    nbrs[2] = ( x_val - 1 + length )%length + y_val * length + z_val * length * length;
    
    /* down, y - 1 */
    nbrs[3] = x_val + (( y_val - 1 + length )%length) * length + z_val * length * length;
    
    /* z + 1 */
    nbrs[4] = x_val  + y_val * length + (( z_val + 1 )%length) * length * length;
    
    /* z - 1 */
    nbrs[5] = x_val  + y_val * length + (( z_val - 1 + length )%length) * length * length;
}


/* PBC in two directions: 
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1).   
 * This results in PBC in the Y and Z directions.
 */
void Cubic::setNbrsYZ(size_t i) 
{

	/* populates the neighbor array */
	size_t x_val = 0, y_val = 0, z_val = 0;
    size_t len2 = length*length;

    /* Get the coordinates of the lattice site */
    x_val = i % length;
    y_val = ( (i - x_val)/length ) % length;
    z_val = (i - x_val - y_val*length) / len2;

    nbrs.clear();

    /* right, x + 1 */
    if( x_val != length - 1 ) {
        nbrs.push_back( ( x_val + 1 ) + y_val * length + z_val * len2 );
    }
    
    /* up, y + 1 */
    nbrs.push_back( x_val + (( y_val + 1 )%length) * length + z_val * len2 );

    /* left, x - 1 */
    if( x_val != 0 ) {
        nbrs.push_back( ( x_val - 1) + y_val * length + z_val * len2 );
    }
    
    /* down, y - 1 */
    nbrs.push_back( x_val + (( y_val - 1 + length )%length) * length + z_val * len2 );
    
    /* z + 1 */
    nbrs.push_back( x_val  + y_val * length + (( z_val + 1 )%length) * len2 );
    
    /* z - 1 */
    nbrs.push_back( x_val  + y_val * length + (( z_val - 1 + length )%length) * len2 );
}


/* PBC in one direction: 
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1).   
 * We do the same for y.
 * This results in PBC in the Z direction.
 */
void Cubic::setNbrsZ(size_t i) 
{

	/* populates the neighbor array */
	size_t x_val = 0, y_val = 0, z_val = 0;
    size_t len2 = length*length;

    /* Get the coordinates of the lattice site */
    x_val = i % length;
    y_val = ( (i - x_val)/length ) % length;
    z_val = (i - x_val - y_val*length) / len2;

    nbrs.clear();

    /* right, x + 1 */
    if( x_val != length - 1 ) {
        nbrs.push_back( ( x_val + 1 ) + y_val * length + z_val * len2 );
    }
    
    /* up, y + 1 */
    if( y_val != length - 1 ) {
        nbrs.push_back( x_val + ( y_val + 1 ) * length + z_val * len2 );
    }

    /* left, x - 1 */
    if( x_val != 0 ) {
        nbrs.push_back( ( x_val - 1) + y_val * length + z_val * len2 );
    }
    
    /* down, y - 1 */
    if( y_val != 0 ) {
        nbrs.push_back( x_val + ( y_val - 1 ) * length + z_val * len2 );
    }
    
    /* z + 1 */
    nbrs.push_back( x_val  + y_val * length + (( z_val + 1 )%length) * len2 );
    
    /* z - 1 */
    nbrs.push_back( x_val  + y_val * length + (( z_val - 1 + length )%length) * len2 );
}


/* PBC in no direction: 
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1).   
 * We do the same for y and z.
 */
void Cubic::setNbrs0(size_t i) 
{

	/* populates the neighbor array */
	size_t x_val = 0, y_val = 0, z_val = 0;
    size_t len2 = length*length;

    /* Get the coordinates of the lattice site */
    x_val = i % length;
    y_val = ( (i - x_val)/length ) % length;
    z_val = (i - x_val - y_val*length) / len2;

    nbrs.clear();

    /* right, x + 1 */
    if( x_val != length - 1 ) {
        nbrs.push_back( ( x_val + 1 ) + y_val * length + z_val * len2 );
    }
    
    /* up, y + 1 */
    if( y_val != length - 1 ) {
        nbrs.push_back( x_val + ( y_val + 1 ) * length + z_val * len2 );
    }

    /* left, x - 1 */
    if( x_val != 0 ) {
        nbrs.push_back( ( x_val - 1) + y_val * length + z_val * len2 );
    }
    
    /* down, y - 1 */
    if( y_val != 0 ) {
        nbrs.push_back( x_val + ( y_val - 1 ) * length + z_val * len2 );
    }
    
    /* z + 1 */
    if( z_val != length - 1 ) {
        nbrs.push_back( x_val  + y_val * length + ( z_val + 1 ) * len2 );
    }
    
    /* z - 1 */
    if( z_val != 0 ) {
        nbrs.push_back( x_val  + y_val * length + ( z_val - 1 ) * len2 );
    }
}
