#include "Cubic4d.h"

/* PBC in all directions */

template<>
size_t Cubic4d<4>::getNumNeighbors(size_t) {
    return 8;
}


template<>
void Cubic4d<4>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t w_val = i % length;
    size_t x_val = ((i - w_val) % length2)/length;
    size_t y_val = ((i - w_val - x_val*length) % length3) / length2;
    size_t z_val = (i - w_val - x_val*length - y_val*length2) / length3;

    /* w + 1 */
    nbrs[0] = (w_val+1)%length + x_val*length + y_val*length2 + z_val*length3;

    /* w - 1 */
    nbrs[1] = (w_val-1+length)%length + x_val*length + y_val*length2 + z_val*length3;

    /* x + 1 */
    nbrs[2] = w_val + ((x_val+1)%length)*length + y_val*length2 + z_val*length3;

    /* x - 1 */
    nbrs[3] = w_val + ((x_val-1+length)%length)*length + y_val*length2 + z_val*length3;

    /* y + 1 */
    nbrs[4] = w_val + x_val*length + ((y_val+1)%length)*length2 + z_val*length3;

    /* y - 1 */
    nbrs[5] = w_val + x_val*length + ((y_val-1+length)%length)*length2 + z_val*length3;

    /* z + 1 */
    nbrs[6] = w_val + x_val*length + y_val*length2 + ((z_val+1)%length)*length3;

    /* z - 1 */
    nbrs[7] = w_val + x_val*length + y_val*length2 + ((z_val-1+length)%length)*length3;
}


/* PBC in three directions: 
 * We check whether the w-value is 0 or length-1. If it is then we don't connect
 * its neighbor across the boundary.
 */
template<>
size_t Cubic4d<3>::getNumNeighbors(size_t i) {
    size_t w_val = i % length;

    size_t nn = 8;
    if(w_val == 0 || w_val == length-1) nn -= 1;
    return nn;
}


template<>
void Cubic4d<3>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t w_val = i % length;
    size_t x_val = ((i - w_val) % length2)/length;
    size_t y_val = ((i - w_val - x_val*length) % length3) / length2;
    size_t z_val = (i - w_val - x_val*length - y_val*length2) / length3;

    if( w_val != length - 1 ) {
        /* w + 1 */
        nbrs[0] = (w_val+1) + x_val*length + y_val*length2 + z_val*length3;
    }

    if( w_val != 0 ) {
        /* w - 1 */
        nbrs[1] = (w_val-1) + x_val*length + y_val*length2 + z_val*length3;
    }

    /* x + 1 */
    nbrs[2] = w_val + ((x_val+1)%length)*length + y_val*length2 + z_val*length3;

    /* x - 1 */
    nbrs[3] = w_val + ((x_val-1+length)%length)*length + y_val*length2 + z_val*length3;

    /* y + 1 */
    nbrs[4] = w_val + x_val*length + ((y_val+1)%length)*length2 + z_val*length3;

    /* y - 1 */
    nbrs[5] = w_val + x_val*length + ((y_val-1+length)%length)*length2 + z_val*length3;

    /* z + 1 */
    nbrs[6] = w_val + x_val*length + y_val*length2 + ((z_val+1)%length)*length3;

    /* z - 1 */
    nbrs[7] = w_val + x_val*length + y_val*length2 + ((z_val-1+length)%length)*length3;
}


/* PBC in two directions: 
 * We check whether the w- or x-value is 0 or length-1. If it is then we don't
 * connect its neighbor across the boundary.
 */
template<>
size_t Cubic4d<2>::getNumNeighbors(size_t i) {
    size_t w_val = i % length;
    size_t x_val = ((i - w_val) % length2)/length;
    
    size_t nn = 8;
    if(w_val == 0 || w_val == length-1) nn -= 1;
    if(x_val == 0 || x_val == length-1) nn -= 1;
    return nn;
}


template<>
void Cubic4d<2>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t w_val = i % length;
    size_t x_val = ((i - w_val) % length2)/length;
    size_t y_val = ((i - w_val - x_val*length) % length3) / length2;
    size_t z_val = (i - w_val - x_val*length - y_val*length2) / length3;

    if( w_val != length - 1 ) {
        /* w + 1 */
        nbrs[0] = (w_val+1) + x_val*length + y_val*length2 + z_val*length3;
    }

    if( w_val != 0 ) {
        /* w - 1 */
        nbrs[1] = (w_val-1) + x_val*length + y_val*length2 + z_val*length3;
    }

    if( x_val != length - 1 ) {
        /* x + 1 */
        nbrs[2] = w_val + (x_val+1)*length + y_val*length2 + z_val*length3;
    }

    if( x_val != 0 ) {
        /* x - 1 */
        nbrs[3] = w_val + (x_val-1)*length + y_val*length2 + z_val*length3;
    }

    /* y + 1 */
    nbrs[4] = w_val + x_val*length + ((y_val+1)%length)*length2 + z_val*length3;

    /* y - 1 */
    nbrs[5] = w_val + x_val*length + ((y_val-1+length)%length)*length2 + z_val*length3;

    /* z + 1 */
    nbrs[6] = w_val + x_val*length + y_val*length2 + ((z_val+1)%length)*length3;

    /* z - 1 */
    nbrs[7] = w_val + x_val*length + y_val*length2 + ((z_val-1+length)%length)*length3;
}

/* PBC in one directions: 
 * We check whether the w-, x-, or y-value is 0 or length-1. If it is then we
 * don't connect its neighbor across the boundary.
 */
template<>
size_t Cubic4d<1>::getNumNeighbors(size_t i) {
    size_t w_val = i % length;
    size_t x_val = ((i - w_val) % length2)/length;
    size_t y_val = ((i - w_val - x_val*length) % length3) / length2;
    
    size_t nn = 8;
    if(w_val == 0 || w_val == length-1) nn -= 1;
    if(x_val == 0 || x_val == length-1) nn -= 1;
    if(y_val == 0 || y_val == length-1) nn -= 1;
    return nn;
}

template<>
void Cubic4d<1>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t w_val = i % length;
    size_t x_val = ((i - w_val) % length2)/length;
    size_t y_val = ((i - w_val - x_val*length) % length3) / length2;
    size_t z_val = (i - w_val - x_val*length - y_val*length2) / length3;

    if( w_val != length - 1 ) {
        /* w + 1 */
        nbrs[0] = (w_val+1) + x_val*length + y_val*length2 + z_val*length3;
    }

    if( w_val != 0 ) {
        /* w - 1 */
        nbrs[1] = (w_val-1) + x_val*length + y_val*length2 + z_val*length3;
    }

    if( x_val != length - 1 ) {
        /* x + 1 */
        nbrs[2] = w_val + (x_val+1)*length + y_val*length2 + z_val*length3;
    }

    if( x_val != 0 ) {
        /* x - 1 */
        nbrs[3] = w_val + (x_val-1)*length + y_val*length2 + z_val*length3;
    }

    if( y_val != length - 1) {
        /* y + 1 */
        nbrs[4] = w_val + x_val*length + (y_val+1)*length2 + z_val*length3;
    }

    if( y_val != 0) {
        /* y - 1 */
        nbrs[5] = w_val + x_val*length + (y_val-1)*length2 + z_val*length3;
    }

    /* z + 1 */
    nbrs[6] = w_val + x_val*length + y_val*length2 + ((z_val+1)%length)*length3;

    /* z - 1 */
    nbrs[7] = w_val + x_val*length + y_val*length2 + ((z_val-1+length)%length)*length3;
}


/* PBC in no directions: 
 * We check whether the w-, x-, y-, or z-value is 0 or length-1. If it is then
 * we don't connect its neighbor across the boundary.
 */
template<>
size_t Cubic4d<0>::getNumNeighbors(size_t i) {
    size_t w_val = i % length;
    size_t x_val = ((i - w_val) % length2)/length;
    size_t y_val = ((i - w_val - x_val*length) % length3) / length2;
    size_t z_val = (i - w_val - x_val*length - y_val*length2) / length3;
    
    size_t nn = 8;
    if(w_val == 0 || w_val == length-1) nn -= 1;
    if(x_val == 0 || x_val == length-1) nn -= 1;
    if(y_val == 0 || y_val == length-1) nn -= 1;
    if(z_val == 0 || z_val == length-1) nn -= 1;
    return nn;
}


template<>
void Cubic4d<0>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t w_val = i % length;
    size_t x_val = ((i - w_val) % length2)/length;
    size_t y_val = ((i - w_val - x_val*length) % length3) / length2;
    size_t z_val = (i - w_val - x_val*length - y_val*length2) / length3;

    if( w_val != length - 1 ) {
        /* w + 1 */
        nbrs[0] = (w_val+1) + x_val*length + y_val*length2 + z_val*length3;
    }

    if( w_val != 0 ) {
        /* w - 1 */
        nbrs[1] = (w_val-1) + x_val*length + y_val*length2 + z_val*length3;
    }

    if( x_val != length - 1 ) {
        /* x + 1 */
        nbrs[2] = w_val + (x_val+1)*length + y_val*length2 + z_val*length3;
    }

    if( x_val != 0 ) {
        /* x - 1 */
        nbrs[3] = w_val + (x_val-1)*length + y_val*length2 + z_val*length3;
    }

    if( y_val != length - 1) {
        /* y + 1 */
        nbrs[4] = w_val + x_val*length + (y_val+1)*length2 + z_val*length3;
    }

    if( y_val != 0) {
        /* y - 1 */
        nbrs[5] = w_val + x_val*length + (y_val-1)*length2 + z_val*length3;
    }

    if( z_val != length - 1) {
        /* z + 1 */
        nbrs[6] = w_val + x_val*length + y_val*length2 + (z_val+1)*length3;
    }

    if( z_val != 0) {
        /* z - 1 */
        nbrs[7] = w_val + x_val*length + y_val*length2 + (z_val-1)*length3;
    }
}


// Let the compiler know we want these
template class Cubic4d<0>;
template class Cubic4d<1>;
template class Cubic4d<2>;
template class Cubic4d<3>;