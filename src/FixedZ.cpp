/* FixedZ.h
 * This class is meant to be a superclass random graphs.  It is meant for site
 * percolation only. A different scheme for keeping track of bonds is necessary
 * for bond-percolation.
 *
 * The neighbor retrieval scheme is a neighbor table.
 *
 * 08/01/05 Created this as a subclass of Lattice.
*/

#include <iostream>

#include "FixedZ.h"
#include "rand.h"
#include "utils.h"

void FixedZ::Setup(size_t z)
{
    max_neighbors = z;
    length = num_sites;

    /* This snippet of code makes the neighbor double array
     * sequential in memory. This speeds things up a bit.
     */
    nbr_ref = new int [num_sites * max_neighbors];
    if( nbr_ref == NULL ) { 
        std::cerr << "Memory alocation error!" << std::endl;
        exit(1);
    }
    for( size_t i=0; i < num_sites * max_neighbors; i++ ) nbr_ref[i] = -1;

    neighbor = new int* [num_sites];
    if( neighbor == NULL ) { 
        std::cerr << "Memory alocation error!" << std::endl;
        exit(1);
    }
    for( size_t j = 0; j < num_sites; j++) {
        neighbor[j] = &nbr_ref[max_neighbors * j];
    }

    // These are needed for bond generation.
    coordination = new size_t [num_sites];
    choose_i = new size_t [num_sites];
    if( coordination == NULL || choose_i == NULL ) { 
        std::cerr << "Memory alocation error!" << std::endl;
        exit(1);
    }
}

void FixedZ::Cleanup()
{
    delete [] neighbor;
    delete [] nbr_ref;
    delete [] coordination;
    delete [] choose_i;
}

size_t FixedZ::getNumNeighbors( size_t )
{
    return max_neighbors;
}


/* When this takes place, a new graph should be generated.
 */
void FixedZ::activateSites() {
    LatticeRandom::activateSites();
    generateBonds();
}

void FixedZ::generateBonds()
{

    size_t i, k, sites_to_choose = num_sites;
    rng_type rng = makeRNG(7);

    for( size_t i = 0; i < num_sites; i++ ) {
        coordination[i] = 0;
        choose_i[i] = i;
    }

    for( size_t bond = 0; bond < max_neighbors; bond++ ) {
        for( size_t count = 0; count < num_sites/2; count++ ) {
            /* pick a random site */
            i = rng(0, sites_to_choose);
            sites_to_choose--;
            /* put that value at the top of the choose list */
            swap( choose_i[i], choose_i[sites_to_choose] );
            i = choose_i[sites_to_choose];

            /* pick a different site */
            k =rng(0, sites_to_choose);
            sites_to_choose--;
            /* put that value at the top of the choose list */
            swap( choose_i[k], choose_i[sites_to_choose] );
            k = choose_i[sites_to_choose];

            neighbor[i][coordination[i]] = k;
            neighbor[k][coordination[k]] = i;
            coordination[i]++;
            coordination[k]++;
        }//count
        sites_to_choose = num_sites;
    }//bond
    //for( size_t i = 0; i < num_sites; i++ ) std::cout << coordination[i] << ' ';

}
