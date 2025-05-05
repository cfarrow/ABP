#include "LatticeRegular.h"
#include <vector>
#include <cmath>

size_t LatticeRegular::getNumNeighbors( size_t i) { 
    if( last != i ) {
        setNbrs(i);
        last = i;
    }
    return nbrs.size();
}

/* For regular graphs, the convention is to create an array of neighbors of site
 * i with this function. It then is called whenever a neighbor is needed. This
 * saves space since a bond array is not required. This will not work with
 * random graphs, however.  THIS IS NOT DEFINED IN LATTICE_H SINCE IT IS
 * SPECIFIC TO REGULAR GRAPHS. See Triangular.h for an example.
 */

/* Returns neighbor number j of site i.  */
size_t LatticeRegular::getNbr(size_t i, size_t j) 
{
    if( last != i ) {
        setNbrs(i);
        last = i;
    }

    return nbrs[j];
}

size_t LatticeRegular::getNumActiveNeighbors(size_t i) 
{

    size_t count = 0;

	if( last != i ) {
        setNbrs(i);
        last = i;
    }

    for( size_t j = 0; j < getNumNeighbors(i); j++ ) {
        if( isActive( nbrs[j] ) ) {
            count++;
        }
    }

    return count;
}