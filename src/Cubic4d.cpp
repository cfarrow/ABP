#include "Cubic4d.h"
#include <vector>
#include <cmath>

void Cubic4d::Setup(short pbc)
{
    /* If num_sites was changed in the constructor, this makes sure that
     * length is properly defined.
     */
    length = static_cast< size_t >( ceil( pow(num_sites, 0.25)) );
    num_neighbors = 8;
    nbrs.resize(num_neighbors);
    last = num_sites;

    /* Set the boundary conditions and set the correct setNbrs function */
    PBC = pbc;
    if( PBC >= 4 ) setNbrs = (LRFptr) &Cubic4d::setNbrsWXYZ;
    else if( PBC == 3 ) setNbrs = (LRFptr) &Cubic4d::setNbrsXYZ;
    else if( PBC == 2 ) setNbrs = (LRFptr) &Cubic4d::setNbrsYZ;
    else if( PBC == 1 ) setNbrs = (LRFptr) &Cubic4d::setNbrsZ;
    else setNbrs = (LRFptr) &Cubic4d::setNbrs0;
}


/* PBC in all directions */
void Cubic4d::setNbrsWXYZ(size_t i) 
{
    /* If periodic boundary conditions are to be added, this is the place! */

	/* populates the neighbor array */
	size_t w_val, x_val, y_val, z_val;
	w_val = x_val = y_val = z_val = 0;

    size_t len2 = length*length;
    size_t len3 = len2*length;

    /* Get the coordinates of the lattice site */
    w_val = i % length;
    x_val = ((i - w_val) % len2)/length;
    y_val = ((i - w_val - x_val*length) % len3) / len2;
    z_val = (i - w_val - x_val*length - y_val*len2) / len3;

    /* w + 1 */
    nbrs[0] = (w_val+1)%length + x_val*length + y_val*len2 + z_val*len3;

    /* w - 1 */
    nbrs[1] = (w_val-1+length)%length + x_val*length + y_val*len2 + z_val*len3;

    /* x + 1 */
    nbrs[2] = w_val + ((x_val+1)%length)*length + y_val*len2 + z_val*len3;

    /* x - 1 */
    nbrs[3] = w_val + ((x_val-1+length)%length)*length + y_val*len2 + z_val*len3;

    /* y + 1 */
    nbrs[4] = w_val + x_val*length + ((y_val+1)%length)*len2 + z_val*len3;

    /* y - 1 */
    nbrs[5] = w_val + x_val*length + ((y_val-1+length)%length)*len2 + z_val*len3;

    /* z + 1 */
    nbrs[6] = w_val + x_val*length + y_val*len2 + ((z_val+1)%length)*len3;

    /* z - 1 */
    nbrs[7] = w_val + x_val*length + y_val*len2 + ((z_val-1+length)%length)*len3;
}


/* PBC in three directions: 
 * We check whether the w-value is 0 or length-1. If it is then we don't connect
 * its neighbor across the boundary.
 */
void Cubic4d::setNbrsXYZ(size_t i) 
{
    /* If periodic boundary conditions are to be added, this is the place! */

	/* populates the neighbor array */
	size_t w_val, x_val, y_val, z_val;
	w_val = x_val = y_val = z_val = 0;

    size_t len2 = length*length;
    size_t len3 = len2*length;

    /* Get the coordinates of the lattice site */
    w_val = i % length;
    x_val = ((i - w_val) % len2)/length;
    y_val = ((i - w_val - x_val*length) % len3) / len2;
    z_val = (i - w_val - x_val*length - y_val*len2) / len3;

    if( w_val != length - 1 ) {
        /* w + 1 */
        nbrs[0] = (w_val+1) + x_val*length + y_val*len2 + z_val*len3;
    }

    if( w_val != 0 ) {
        /* w - 1 */
        nbrs[1] = (w_val-1) + x_val*length + y_val*len2 + z_val*len3;
    }

    /* x + 1 */
    nbrs[2] = w_val + ((x_val+1)%length)*length + y_val*len2 + z_val*len3;

    /* x - 1 */
    nbrs[3] = w_val + ((x_val-1+length)%length)*length + y_val*len2 + z_val*len3;

    /* y + 1 */
    nbrs[4] = w_val + x_val*length + ((y_val+1)%length)*len2 + z_val*len3;

    /* y - 1 */
    nbrs[5] = w_val + x_val*length + ((y_val-1+length)%length)*len2 + z_val*len3;

    /* z + 1 */
    nbrs[6] = w_val + x_val*length + y_val*len2 + ((z_val+1)%length)*len3;

    /* z - 1 */
    nbrs[7] = w_val + x_val*length + y_val*len2 + ((z_val-1+length)%length)*len3;
}


/* PBC in two directions: 
 * We check whether the w- or x-value is 0 or length-1. If it is then we don't
 * connect its neighbor across the boundary.
 */
void Cubic4d::setNbrsYZ(size_t i) 
{
    /* If periodic boundary conditions are to be added, this is the place! */

	/* populates the neighbor array */
	size_t w_val, x_val, y_val, z_val;
	w_val = x_val = y_val = z_val = 0;

    size_t len2 = length*length;
    size_t len3 = len2*length;

    /* Get the coordinates of the lattice site */
    w_val = i % length;
    x_val = ((i - w_val) % len2)/length;
    y_val = ((i - w_val - x_val*length) % len3) / len2;
    z_val = (i - w_val - x_val*length - y_val*len2) / len3;

    if( w_val != length - 1 ) {
        /* w + 1 */
        nbrs[0] = (w_val+1) + x_val*length + y_val*len2 + z_val*len3;
    }

    if( w_val != 0 ) {
        /* w - 1 */
        nbrs[1] = (w_val-1) + x_val*length + y_val*len2 + z_val*len3;
    }

    if( x_val != length - 1 ) {
        /* x + 1 */
        nbrs[2] = w_val + (x_val+1)*length + y_val*len2 + z_val*len3;
    }

    if( x_val != 0 ) {
        /* x - 1 */
        nbrs[3] = w_val + (x_val-1)*length + y_val*len2 + z_val*len3;
    }

    /* y + 1 */
    nbrs[4] = w_val + x_val*length + ((y_val+1)%length)*len2 + z_val*len3;

    /* y - 1 */
    nbrs[5] = w_val + x_val*length + ((y_val-1+length)%length)*len2 + z_val*len3;

    /* z + 1 */
    nbrs[6] = w_val + x_val*length + y_val*len2 + ((z_val+1)%length)*len3;

    /* z - 1 */
    nbrs[7] = w_val + x_val*length + y_val*len2 + ((z_val-1+length)%length)*len3;
}

/* PBC in one directions: 
 * We check whether the w-, x-, or y-value is 0 or length-1. If it is then we
 * don't connect its neighbor across the boundary.
 */
void Cubic4d::setNbrsZ(size_t i) 
{
    /* If periodic boundary conditions are to be added, this is the place! */

	/* populates the neighbor array */
	size_t w_val, x_val, y_val, z_val;
	w_val = x_val = y_val = z_val = 0;

    size_t len2 = length*length;
    size_t len3 = len2*length;

    /* Get the coordinates of the lattice site */
    w_val = i % length;
    x_val = ((i - w_val) % len2)/length;
    y_val = ((i - w_val - x_val*length) % len3) / len2;
    z_val = (i - w_val - x_val*length - y_val*len2) / len3;

    if( w_val != length - 1 ) {
        /* w + 1 */
        nbrs[0] = (w_val+1) + x_val*length + y_val*len2 + z_val*len3;
    }

    if( w_val != 0 ) {
        /* w - 1 */
        nbrs[1] = (w_val-1) + x_val*length + y_val*len2 + z_val*len3;
    }

    if( x_val != length - 1 ) {
        /* x + 1 */
        nbrs[2] = w_val + (x_val+1)*length + y_val*len2 + z_val*len3;
    }

    if( x_val != 0 ) {
        /* x - 1 */
        nbrs[3] = w_val + (x_val-1)*length + y_val*len2 + z_val*len3;
    }

    if( y_val != length - 1) {
        /* y + 1 */
        nbrs[4] = w_val + x_val*length + (y_val+1)*len2 + z_val*len3;
    }

    if( y_val != 0) {
        /* y - 1 */
        nbrs[5] = w_val + x_val*length + (y_val-1)*len2 + z_val*len3;
    }

    /* z + 1 */
    nbrs[6] = w_val + x_val*length + y_val*len2 + ((z_val+1)%length)*len3;

    /* z - 1 */
    nbrs[7] = w_val + x_val*length + y_val*len2 + ((z_val-1+length)%length)*len3;
}


/* PBC in no directions: 
 * We check whether the w-, x-, y-, or z-value is 0 or length-1. If it is then
 * we don't connect its neighbor across the boundary.
 */
void Cubic4d::setNbrs0(size_t i) 
{
    /* If periodic boundary conditions are to be added, this is the place! */

	/* populates the neighbor array */
	size_t w_val, x_val, y_val, z_val;
	w_val = x_val = y_val = z_val = 0;

    size_t len2 = length*length;
    size_t len3 = len2*length;

    /* Get the coordinates of the lattice site */
    w_val = i % length;
    x_val = ((i - w_val) % len2)/length;
    y_val = ((i - w_val - x_val*length) % len3) / len2;
    z_val = (i - w_val - x_val*length - y_val*len2) / len3;

    if( w_val != length - 1 ) {
        /* w + 1 */
        nbrs[0] = (w_val+1) + x_val*length + y_val*len2 + z_val*len3;
    }

    if( w_val != 0 ) {
        /* w - 1 */
        nbrs[1] = (w_val-1) + x_val*length + y_val*len2 + z_val*len3;
    }

    if( x_val != length - 1 ) {
        /* x + 1 */
        nbrs[2] = w_val + (x_val+1)*length + y_val*len2 + z_val*len3;
    }

    if( x_val != 0 ) {
        /* x - 1 */
        nbrs[3] = w_val + (x_val-1)*length + y_val*len2 + z_val*len3;
    }

    if( y_val != length - 1) {
        /* y + 1 */
        nbrs[4] = w_val + x_val*length + (y_val+1)*len2 + z_val*len3;
    }

    if( y_val != 0) {
        /* y - 1 */
        nbrs[5] = w_val + x_val*length + (y_val-1)*len2 + z_val*len3;
    }

    if( z_val != length - 1) {
        /* z + 1 */
        nbrs[6] = w_val + x_val*length + y_val*len2 + (z_val+1)*len3;
    }

    if( z_val != 0) {
        /* z - 1 */
        nbrs[7] = w_val + x_val*length + y_val*len2 + (z_val-1)*len3;
    }
}
