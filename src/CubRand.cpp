#include <iostream>

#include "CubRand.h"
#include "rand.h"
#include "utils.h"

void CubRand::Setup(double fb)
{

    if(fb > 1) fb = 1;
    if(fb < 0) fb = 0;
    f = fb;
    nb = static_cast<size_t>(3*f*num_sites);
    nb -= nb%2; // Make sure this is even!

    max_neighbors = 6;
    /* If num_sites was changed in the constructor, this makes sure that
     * length is properly defined.
     */
    length = static_cast< size_t >( ceil( pow(num_sites, 1.0/3.0) ) );

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
    bondlist = new int [3*num_sites];
    if( bondlist == NULL) {
        std::cerr << "Memory alocation error!" << std::endl;
        exit(1);
    }
    for(size_t l=0; l < 3*num_sites; ++l) {
        bondlist[l] = l;
    }

}

void CubRand::Cleanup()
{
    delete [] neighbor;
    delete [] nbr_ref;
    delete [] bondlist;
}

/* When this takes place, a new graph should be generated.
 */
void CubRand::activateSites() {
    LatticeRandom::activateSites();
    generateBonds();
}

void CubRand::generateBonds() {

    // Start by making a cubic lattice
    size_t i, k;
	size_t x_val, y_val, z_val;

    for(i=0; i<num_sites; ++i) {
        /* Get the coordinates of the lattice site */
        x_val = i % length;
        y_val = ( (i - x_val)/length ) % length;
        z_val = (i - x_val - y_val*length) / ( length * length );

        /* right, x + 1 */
        k = ( x_val + 1 )%length + y_val * length + z_val * length * length;
        neighbor[i][0] = k;
        neighbor[k][3] = i;
        
        /* up, y + 1 */
        k = x_val + (( y_val + 1 )%length) * length + z_val * length * length;
        neighbor[i][1] = k;
        neighbor[k][4] = i;
        
        /* z + 1 */
        k = x_val  + y_val * length + (( z_val + 1 )%length) * length * length;
        neighbor[i][2] = k;
        neighbor[k][5] = i;
    }

    /* Remove a fraction f of bonds (nb bonds total).
     *
     * Note that any bond can be indexed as follows:
     * nbr = index % 3
     * site = (index - nbr)/3
     *
     * Thus, the head and tail of the bond are:
     * site, neighbor[site][nbr]
     */

    // Randomly shuffle the bonds to be removed to the back of the bondlist. The
    // last nb bonds will be those that need to be replaced. Note that there is
    // no possiblity of removing the same bond twice with this scheme. In
    // addition, it can efficiently remove a large number of bonds. The tradeoff
    // is the necessity of a bond list.
    //
    rng_type rng = makeRNG(7);
    int bidx1;
    size_t bonds_left = 3*num_sites;

    for(size_t l=0; l < nb; ++l) {

        bidx1 = rng(0, bonds_left);
        bonds_left--;
        swap( bondlist[bidx1], bondlist[bonds_left] );
    }

    /* Now randomly reconnect the bonds in the list. Do this by cycling through
     * the removed bonds. Grab two random bonds and connect them at the tails.
     * After cycling through the bonds once, do it again and connect random
     * pairs at the heads.
     * 
     * At each step, it must be checked that no self-loops and no double-bonds
     * are created.
     */

    int nidx1, nidx2, bidx2;
    size_t bonds_to_choose;
    bool step_back;

    for(size_t loop = 0; loop < 2; loop++)
    {
        bonds_to_choose = nb;
        for(int l=0; l < nb/2; ++l)
        {

            bidx1 = rng(0, bonds_to_choose);
            bidx1 += bonds_left;
            bonds_to_choose--;
            swap( bondlist[bidx1], bondlist[bonds_left + bonds_to_choose] );
            bidx1 = bondlist[bonds_left + bonds_to_choose];

            bidx2 = rng(0, bonds_to_choose);
            bidx2 += bonds_left;
            bonds_to_choose--;
            swap( bondlist[bidx2], bondlist[bonds_left + bonds_to_choose] );
            bidx2 = bondlist[bonds_left + bonds_to_choose];

            /* connect the bonds */

            // site i
            nidx1 = bidx1 % 3;
            i = (bidx1 - nidx1)/3;
            if( loop == 0 )
            {
                i = neighbor[i][nidx1];
                nidx1 = (nidx1 + 3) % 6;
            }

            // site k
            nidx2 = bidx2 % 3;
            k = (bidx2 - nidx2)/3;
            if( loop == 0 )
            {
                k = neighbor[k][nidx2];
                nidx2 = (nidx2 + 3) % 6;
            }

            // Check for self-loops and double-bonds.
            // This is a bad hack and needs to be improved.
            step_back = 0;
            if(i==k) step_back = 1;
            for(int j=0; j<6; ++j) {
                if(j != nidx1 && neighbor[i][j] == k) {
                    step_back = 1;
                    break;
                }
            }
            if(step_back) {
                bonds_to_choose += 2;
                --l;
                // Start over if we're going to get stuck.
                if(bonds_to_choose <= 6) {
                    bonds_to_choose = nb;
                    l = -1;
                }
                continue;
            }

            // Connect the sites.
            neighbor[i][nidx1] = k;
            neighbor[k][nidx2] = i;
        }
    }
}
