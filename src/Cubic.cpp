#include "Cubic.h"

/* PBC in all directions */

template<>
size_t Cubic<3>::getNumNeighbors(size_t) {
    return 6;
}


template <>
void Cubic<3>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( (i - x_val)/length ) % length;
    size_t z_val = (i - x_val - y_val*length) / ( length * length );

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
template<>
size_t Cubic<2>::getNumNeighbors(size_t i) {
    size_t x_val = i % length;
    size_t nn = 6;
    if(x_val == 0 || x_val == length-1) nn -= 1;
    return nn;
}


template <>
void Cubic<2>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( (i - x_val)/length ) % length;
    size_t z_val = (i - x_val - y_val*length) / length2;

    nbrs.clear();

    /* right, x + 1 */
    if( x_val != length - 1 ) {
        nbrs.push_back( ( x_val + 1 ) + y_val * length + z_val * length2 );
    }
    
    /* up, y + 1 */
    nbrs.push_back( x_val + (( y_val + 1 )%length) * length + z_val * length2 );

    /* left, x - 1 */
    if( x_val != 0 ) {
        nbrs.push_back( ( x_val - 1) + y_val * length + z_val * length2 );
    }
    
    /* down, y - 1 */
    nbrs.push_back( x_val + (( y_val - 1 + length )%length) * length + z_val * length2 );
    
    /* z + 1 */
    nbrs.push_back( x_val  + y_val * length + (( z_val + 1 )%length) * length2 );
    
    /* z - 1 */
    nbrs.push_back( x_val  + y_val * length + (( z_val - 1 + length )%length) * length2 );
}


/* PBC in one direction: 
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1).   
 * We do the same for y.
 * This results in PBC in the Z direction.
 */
template<>
size_t Cubic<1>::getNumNeighbors(size_t i) {
    size_t x_val = i % length;
    size_t y_val = ( (i - x_val)/length ) % length;
    size_t nn = 6;
    if(x_val == 0 || x_val == length-1) nn -= 1;
    if(y_val == 0 || y_val == length-1) nn -= 1;
    return nn;
}


template <>
void Cubic<1>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( (i - x_val)/length ) % length;
    size_t z_val = (i - x_val - y_val*length) / length2;

    nbrs.clear();

    /* right, x + 1 */
    if( x_val != length - 1 ) {
        nbrs.push_back( ( x_val + 1 ) + y_val * length + z_val * length2 );
    }
    
    /* up, y + 1 */
    if( y_val != length - 1 ) {
        nbrs.push_back( x_val + ( y_val + 1 ) * length + z_val * length2 );
    }

    /* left, x - 1 */
    if( x_val != 0 ) {
        nbrs.push_back( ( x_val - 1) + y_val * length + z_val * length2 );
    }
    
    /* down, y - 1 */
    if( y_val != 0 ) {
        nbrs.push_back( x_val + ( y_val - 1 ) * length + z_val * length2 );
    }
    
    /* z + 1 */
    nbrs.push_back( x_val  + y_val * length + (( z_val + 1 )%length) * length2 );
    
    /* z - 1 */
    nbrs.push_back( x_val  + y_val * length + (( z_val - 1 + length )%length) * length2 );
}


/* PBC in no direction: 
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1).   
 * We do the same for y and z.
 */
template<>
size_t Cubic<0>::getNumNeighbors(size_t i) {
    size_t x_val = i % length;
    size_t y_val = ( (i - x_val)/length ) % length;
    size_t z_val = (i - x_val - y_val*length) / length2;
    size_t nn = 6;
    if(x_val == 0 || x_val == length-1) nn -= 1;
    if(y_val == 0 || y_val == length-1) nn -= 1;
    if(z_val == 0 || z_val == length-1) nn -= 1;
    return nn;
}


template <>
void Cubic<0>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( (i - x_val)/length ) % length;
    size_t z_val = (i - x_val - y_val*length) / length2;

    nbrs.clear();

    /* right, x + 1 */
    if( x_val != length - 1 ) {
        nbrs.push_back( ( x_val + 1 ) + y_val * length + z_val * length2 );
    }
    
    /* up, y + 1 */
    if( y_val != length - 1 ) {
        nbrs.push_back( x_val + ( y_val + 1 ) * length + z_val * length2 );
    }

    /* left, x - 1 */
    if( x_val != 0 ) {
        nbrs.push_back( ( x_val - 1) + y_val * length + z_val * length2 );
    }
    
    /* down, y - 1 */
    if( y_val != 0 ) {
        nbrs.push_back( x_val + ( y_val - 1 ) * length + z_val * length2 );
    }
    
    /* z + 1 */
    if( z_val != length - 1 ) {
        nbrs.push_back( x_val  + y_val * length + ( z_val + 1 ) * length2 );
    }
    
    /* z - 1 */
    if( z_val != 0 ) {
        nbrs.push_back( x_val  + y_val * length + ( z_val - 1 ) * length2 );
    }
}


// Let the compiler know we want these
template class Cubic<0>;
template class Cubic<1>;
template class Cubic<2>;
template class Cubic<3>;