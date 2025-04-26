/* BCC.h
 * 03/06/07 Created this as a subclass of LatticeRegular3d.
*/

#include "BCC.h"
#include <vector>
#include <cmath>

void BCC::Setup( short pbc )
{
    /* If num_sites was changed in the constructor, this makes sure that
     * length is properly defined.
     */
    length = static_cast< size_t >( ceil( pow(num_sites, 1.0/3.0)) );
    num_neighbors = 8;
    nbrs.resize(num_neighbors);
    last = num_sites;

    /* Set the boundary conditions and set the correct setNbrs function */
    PBC = pbc;
    if( PBC >= 3 ) setNbrs = (LRFptr) &BCC::setNbrsXYZ;
    else if( PBC == 2 ) setNbrs = (LRFptr) &BCC::setNbrsYZ;
    else if( PBC == 1 ) setNbrs = (LRFptr) &BCC::setNbrsZ;
    else setNbrs = (LRFptr) &BCC::setNbrs0;

}


/* PBC in all directions */
void BCC::setNbrsXYZ(size_t i) 
{
    /* If periodic boundary conditions are to be added, this is the place! */

	/* populates the neighbor array */
	size_t x_val, y_val, z_val;
    static size_t len2 = length*length;
	x_val = y_val = z_val = 0;

    /* Get the coordinates of the lattice site */
    x_val = i % length;
    y_val = ( (i - x_val)/length ) % length;
    z_val = (i - x_val - y_val*length) / ( length * length );

    /* Consider the inner bcc sublattice placed on the same half-spaced grid as
     * the other. The bcc bonds are not all equal length, but this does not
     * matter in percolation.
     *
     * If on sublattice 'A' (z odd) the neighbors are:
     * x        y       z+/-1
     * x        y+1     z+/-1
     * x+1      y       z+/-1
     * x+1      y+1     z+/-1
     *
     * If on sublattice 'B' (z even) the neighbors are:
     * x-1      y-1     z+/-1
     * x-1      y       z+/-1
     * x        y-1     z+/-1
     * x        y       z+/-1
     */

    /* x, y, z+1 */
    nbrs[0] = x_val + y_val*length + ((z_val+1)%length)*len2;
    
    /* x, y, z-1 */
    nbrs[1] = x_val + y_val*length + ((z_val-1+length)%length)*len2;

    if( z_val %2 )
    {
        /* x, y+1, z+1 */
        nbrs[2] = x_val + ((y_val+1)%length)*length + ((z_val+1)%length)*len2;
        
        /* x, y+1, z-1 */
        nbrs[3] = x_val + ((y_val+1)%length)*length + ((z_val-1+length)%length)*len2;

        /* x+1, y, z+1 */
        nbrs[4] = (x_val+1)%length + y_val*length + ((z_val+1)%length)*len2;
        
        /* x+1, y, z-1 */
        nbrs[5] = (x_val+1)%length + y_val*length + ((z_val-1+length)%length)*len2;
        
        /* x+1, y+1, z+1 */
        nbrs[6] = (x_val+1)%length + ((y_val+1)%length)*length + ((z_val+1)%length)*len2;
        
        /* x+1, y+1, z-1 */
        nbrs[7] = (x_val+1)%length + ((y_val+1)%length)*length + ((z_val-1+length)%length)*len2;

    }

    else
    {
        /* x-1, y-1, z+1 */
        nbrs[2] = (x_val-1+length)%length + ((y_val-1+length)%length)*length + ((z_val+1)%length)*len2;

        /* x-1, y-1, z-1 */
        nbrs[3] = (x_val-1+length)%length + ((y_val-1+length)%length)*length + ((z_val-1+length)%length)*len2;
        
        /* x-1, y, z+1 */
        nbrs[4] = (x_val-1+length)%length + y_val*length + ((z_val+1)%length)*len2;
        
        /* x-1, y, z-1 */
        nbrs[5] = (x_val-1+length)%length + y_val*length + ((z_val-1+length)%length)*len2;

        /* x, y-1, z+1 */
        nbrs[6] = x_val + ((y_val-1+length)%length)*length + ((z_val+1)%length)*len2;
        
        /* x, y-1, z-1 */
        nbrs[7] = x_val + ((y_val-1+length)%length)*length + ((z_val-1+length)%length)*len2;

    }
    
}

/* PBC in two directions: 
 * If x_val == 0, don't connect to "x_val-1" neighbors.
 * If x_val == length - 1, don't connect to "x_val+1" neighbors.
 */
void BCC::setNbrsYZ(size_t i) 
{

	/* populates the neighbor array */
	size_t x_val = 0, y_val = 0, z_val = 0;
    static size_t len2 = length*length;

    /* Get the coordinates of the lattice site */
    x_val = i % length;
    y_val = ( (i - x_val)/length ) % length;
    z_val = (i - x_val - y_val*length) / len2;

    nbrs.clear();

    /* x, y, z+1 */
    nbrs.push_back( x_val + y_val*length + ((z_val+1)%length)*len2 );
    
    /* x, y, z-1 */
    nbrs.push_back( x_val + y_val*length + ((z_val-1+length)%length)*len2 );

    if( z_val %2 )
    {
        /* x, y+1, z+1 */
        nbrs.push_back( x_val + ((y_val+1)%length)*length + ((z_val+1)%length)*len2 );
        
        /* x, y+1, z-1 */
        nbrs.push_back( x_val + ((y_val+1)%length)*length + ((z_val-1+length)%length)*len2 );

        if( x_val != length - 1)
        {
            /* x+1, y, z+1 */
            nbrs.push_back( x_val + y_val*length + ((z_val+1)%length)*len2 );
            
            /* x+1, y, z-1 */
            nbrs.push_back( x_val + y_val*length + ((z_val-1+length)%length)*len2 );
            
            /* x+1, y+1, z+1 */
            nbrs.push_back( x_val + ((y_val+1)%length)*length + ((z_val+1)%length)*len2 );
            
            /* x+1, y+1, z-1 */
            nbrs.push_back( x_val + ((y_val+1)%length)*length + ((z_val-1+length)%length)*len2 );
        }

    }

    else
    {

        if( x_val != 0 )
        {
            /* x-1, y-1, z+1 */
            nbrs.push_back( x_val-1 + ((y_val-1+length)%length)*length + ((z_val+1)%length)*len2 );

            /* x-1, y-1, z-1 */
            nbrs.push_back( x_val-1 + ((y_val-1+length)%length)*length + ((z_val-1+length)%length)*len2 );
            
            /* x-1, y, z+1 */
            nbrs.push_back( x_val-1 + y_val*length + ((z_val+1)%length)*len2 );
            
            /* x-1, y, z-1 */
            nbrs.push_back( x_val-1 + y_val*length + ((z_val-1+length)%length)*len2 );
        }

        /* x, y-1, z+1 */
        nbrs.push_back( x_val + ((y_val-1+length)%length)*length + ((z_val+1)%length)*len2 );
        
        /* x, y-1, z-1 */
        nbrs.push_back( x_val + ((y_val-1+length)%length)*length + ((z_val-1+length)%length)*len2 );

    }

}

/* PBC in one direction: 
 * PBC in two directions: 
 * If x_val == 0, don't connect to "x_val-1" neighbors.
 * If x_val == length - 1, don't connect to "x_val+1" neighbors.
 * If y_val == 0, don't connect to "y_val-1" neighbors.
 * If y_val == length - 1, don't connect to "y_val+1" neighbors.
 */
void BCC::setNbrsZ(size_t i) 
{

	/* populates the neighbor array */
	size_t x_val = 0, y_val = 0, z_val = 0;
    static size_t len2 = length*length;

    /* Get the coordinates of the lattice site */
    x_val = i % length;
    y_val = ( (i - x_val)/length ) % length;
    z_val = (i - x_val - y_val*length) / len2;

    nbrs.clear();

    /* x, y, z+1 */
    nbrs.push_back( x_val + y_val*length + ((z_val+1)%length)*len2 );
    
    /* x, y, z-1 */
    nbrs.push_back( x_val + y_val*length + ((z_val-1+length)%length)*len2 );

    if( z_val %2 )
    {
        if( y_val != length - 1 )
        {
            /* x, y+1, z+1 */
            nbrs.push_back( x_val + (y_val+1)*length + ((z_val+1)%length)*len2 );
            
            /* x, y+1, z-1 */
            nbrs.push_back( x_val + (y_val+1)*length + ((z_val-1+length)%length)*len2 );
        }

        if( x_val != length - 1)
        {
            /* x+1, y, z+1 */
            nbrs.push_back( x_val + y_val*length + ((z_val+1)%length)*len2 );
            
            /* x+1, y, z-1 */
            nbrs.push_back( x_val + y_val*length + ((z_val-1+length)%length)*len2 );
            
            if( y_val != length - 1 )
            {
                /* x+1, y+1, z+1 */
                nbrs.push_back( x_val + (y_val+1)*length + ((z_val+1)%length)*len2 );
                
                /* x+1, y+1, z-1 */
                nbrs.push_back( x_val + (y_val+1)*length + ((z_val-1+length)%length)*len2 );
            }
        }

    }

    else
    {

        if( x_val != 0 )
        {
            if( y_val != 0 ) {
                /* x-1, y-1, z+1 */
                nbrs.push_back( x_val-1 + (y_val-1)*length + ((z_val+1)%length)*len2 );

                /* x-1, y-1, z-1 */
                nbrs.push_back( x_val-1 + (y_val-1)*length + ((z_val-1+length)%length)*len2 );
                
            }

            /* x-1, y, z+1 */
            nbrs.push_back( x_val-1 + y_val*length + ((z_val+1)%length)*len2 );
            
            /* x-1, y, z-1 */
            nbrs.push_back( x_val-1 + y_val*length + ((z_val-1+length)%length)*len2 );
        }

        if( y_val != 0 ) {
            /* x, y-1, z+1 */
            nbrs.push_back( x_val + (y_val-1)*length + ((z_val+1)%length)*len2 );
            
            /* x, y-1, z-1 */
            nbrs.push_back( x_val + (y_val-1)*length + ((z_val-1+length)%length)*len2 );

        }
    }
}


/* PBC in no direction: 
 * We check whether the x-value is 0 or length-1. If it is then we don't
 * connect the neighbor to the left (for 0) or right (for length-1).   
 * We do the same for y and z.
 */
void BCC::setNbrs0(size_t i) 
{

	/* populates the neighbor array */
	size_t x_val = 0, y_val = 0, z_val = 0;
    static size_t len2 = length*length;

    /* Get the coordinates of the lattice site */
    x_val = i % length;
    y_val = ( (i - x_val)/length ) % length;
    z_val = (i - x_val - y_val*length) / len2;

    nbrs.clear();

    if( z_val != length + 1 )
    {
        /* x, y, z+1 */
        nbrs.push_back( x_val + y_val*length + (z_val+1)*len2 );
    }
    
    if( z_val != 0 )
    {
        /* x, y, z-1 */
        nbrs.push_back( x_val + y_val*length + (z_val-1)*len2 );
    }

    if( z_val %2 )
    {
        if(y_val != length-1 && z_val != length -1)
        {
            /* x, y+1, z+1 */
            nbrs.push_back( x_val + (y_val+1)*length + (z_val+1)*len2 );
        }
        
        if(y_val != length-1 && z_val != 0)
        {
            /* x, y+1, z-1 */
            nbrs.push_back( x_val + (y_val+1)*length + (z_val-1)*len2 );
        }

        if(x_val != length-1 && z_val != length - 1)
        {
            /* x+1, y, z+1 */
            nbrs.push_back( x_val+1 + y_val*length + (z_val+1)*len2 );
        }
        
        if(x_val != length-1 && z_val != 0)
        {
            /* x+1, y, z-1 */
            nbrs.push_back( x_val+1 + y_val*length + (z_val-1)*len2 );
        }
        
        if(x_val != length-1 && y_val != length-1 && z_val != length-1)
        {
            /* x+1, y+1, z+1 */
            nbrs.push_back( x_val+1 + (y_val+1)*length + (z_val+1)*len2 );
        }
        
        if(x_val != length-1 && y_val != length-1 && z_val != 0)
        {
            /* x+1, y+1, z-1 */
            nbrs.push_back( x_val+1 + (y_val+1)*length + (z_val-1)*len2 );
        }

    }

    else
    {
        if(x_val != 0 && y_val != 0 && z_val != length-1)
        {
            /* x-1, y-1, z+1 */
            nbrs.push_back( x_val-1 + (y_val-1)*length + (z_val+1)*len2 );
        }

        if(x_val != 0 && y_val != 0 && z_val != 0)
        {
            /* x-1, y-1, z-1 */
            nbrs.push_back( x_val-1 + (y_val-1)*length + (z_val-1)*len2 );
        }
        
        if(x_val != 0 && z_val != length-1)
        {
            /* x-1, y, z+1 */
            nbrs.push_back( x_val-1 + y_val*length + (z_val+1)*len2 );
        }
        
        if(x_val != 0 && z_val != 0)
        {
            /* x-1, y, z-1 */
            nbrs.push_back( x_val-1 + y_val*length + (z_val-1)*len2 );
        }

        if(y_val != 0 && z_val != length-1)
        {
            /* x, y-1, z+1 */
            nbrs.push_back( x_val + (y_val-1)*length + (z_val+1)*len2 );
        }
        
        if(y_val != 0 && z_val != 0)
        {
            /* x, y-1, z-1 */
            nbrs.push_back( x_val + (y_val-1)*length + (z_val-1)*len2 );
        }

    }
    
}
