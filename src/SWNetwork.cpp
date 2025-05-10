#include <iostream>
#include <random>

#include "SWNetwork.h"
#include "rand.h"
#include "utils.h"

// TODO

void SWNetwork::Setup(size_t _dims, double _alpha, double _gamma)
{
    // constants
    max_neighbors = 50;
    SHORT_BOND_LENGTH = 2.3;

    dims = _dims;
    alpha = _alpha;
    gamma = _gamma;

    // This will not have roundoff error. See the constructor for details.
    length = static_cast<size_t>(pow(num_sites, 1.0/dims) );

    /* This snippet of code makes the neighbor double array
     * sequential in memory. This speeds things up a bit.
     */
    nbr_ref = new size_t [num_sites * max_neighbors];
    if( nbr_ref == NULL ) { 
        std::cerr << "Memory alocation error!" << std::endl;
        exit(1);
    }
    for( size_t i=0; i < num_sites * max_neighbors; i++ ) 
    {
        nbr_ref[i] = -1;
    }

    neighbor = new size_t* [num_sites];
    if( neighbor == NULL ) { 
        std::cerr << "Memory alocation error!" << std::endl;
        exit(1);
    }
    for( size_t j = 0; j < num_sites; j++) {
        neighbor[j] = &nbr_ref[max_neighbors * j];
    }

    coordination = new size_t [num_sites];
    if( coordination == NULL ) { 
        std::cerr << "Memory alocation error!" << std::endl;
        exit(1);
    }

    bond_weight = new double [num_sites];
    if ( bond_weight == NULL ) {
        std::cerr << "Memory alocation error!" << std::endl;
        exit(1);
    }

    /* Calculate "coordinate jumps" for labelling of regular d-dim. lattice with
     * helical BC */
    jump = new size_t[dims];
    if ( jump == NULL ) {
        std::cerr << "Memory alocation error!" << std::endl;
        exit(1);
    }
    jump[0] = 1;
    for(size_t dir=1; dir < dims; dir++){
        jump[dir]=length*jump[dir-1];
    }

    // calculate the bond weights
    calculateBondWeights();
}

void SWNetwork::Cleanup()
{
    delete [] neighbor;
    delete [] nbr_ref;
    delete [] coordination;
    delete [] bond_weight;
    delete [] jump;
}

/* When this takes place, a new graph should be generated.
 */
void SWNetwork::activateSites() {
    LatticeRandom::activateSites();
    generateBonds();
}

// Code by Cristian Moukarzel, adopted by Chris Farrow
void SWNetwork::generateBonds() {

    size_t inbr;
    size_t isite,ksite,b_label;
    size_t dir,kounter;
    size_t bond;

    size_t half_length = length/2;

	static std::mt19937_64 gen(7);
	static std::uniform_real_distribution<double> dist(0.0, 1.0);
	auto prng = [&]() { return dist(gen); };

    bounded_rng_type rng = makeRNG(7, 0, num_sites-1);

    // Separation of bonds into short and long ones.
    double bond_len, bond_mod;
    std::vector<double> vctr_coord(dims+1,0.0);

    double total_weight,psum,rweight;

    // Parameters for creation of the bond length probability distribution
    // obtained via inversion method.  half_length appears because of
    // the distribution being finite.  Else, the integrals would not converge in
    // case of exponent being >= (-1).
    double ddim = static_cast<double>(dims);
    double dexponent = 1.0/(ddim - alpha);
    double lambda = pow((1.0*half_length + 1.0), ddim - alpha) - 1.0;
    
    // Initialize the coordination array
    for(size_t i = 0; i < num_sites; ++i) {
        coordination[i] = 0;
    }

    // The AVERAGE coordination is 2*bonds per site. Individual sites have a
    // Poisson-distributed number of links (for large systems).

    size_t total_bonds= static_cast<size_t>(0.5 * gamma * num_sites);

    // Now connect a total of total_bonds bonds
    for(bond=0;bond<total_bonds;bond++) { 

        kounter = 0;
        do{ // this loop is repeated until a pair of neighbors is found that
            // can be connected.
            kounter++;

            // NOTE: A new starting site is chosen at each trial 
            isite = rng();

            // Now choose the OTHER site.
            if( prng() < p_s) {   
                // use method for short bonds Loop over short bonds and pick one
                // with probability short_bond_weight[]
                total_weight = 0.0;
                for(b_label=0; b_label < short_bond_label.size(); b_label++){
                    ksite=(isite+short_bond_label[b_label]+num_sites)%num_sites;
                    // Ignore if this site is already connected to isite. 
                    for(inbr=0; inbr < coordination[isite]; inbr++) {
                        if(neighbor[isite][inbr]==ksite) {
                            break;
                        }
                    }
                    if(inbr < coordination[isite]) {
                        // sites are connected. Ignore
                        short_bond_weight[b_label]=0.0;
                    }
                    else { 
                        // site is elligible. 
                        short_bond_weight[b_label] = bond_weight[short_bond_label[b_label]];
                    }
                    total_weight += short_bond_weight[b_label];
                } // end loop over short bonds

                if(total_weight > 0.0) {
                    // Normalize and choose one
                    rweight = prng();
                    psum = 0.0;
                    for(b_label=0; b_label < short_bond_label.size(); b_label++) {
                        short_bond_weight[b_label] /= total_weight;
                        psum += short_bond_weight[b_label];
                        if(rweight < psum) {
                            ksite=(isite+short_bond_label[b_label]+num_sites)%num_sites;
                            break;
                        }
                    }
                    break;
                }
                else {
                    // no linkable site found
                    continue; // jump to end of loop - a new isite will be chosen
                }	
            }
            else { // use method for long bonds
                // The LONG-BOND method consists of the folowing steps: 1) picking
                // a random distance, 2) picking a random dir, and then
                // choosing the site closest to the point so generated.
                
                // 1) Generate bond_len, a real random variable between 1 and
                // L/2, with PowerLaw distribution.
                do {
                    if(alpha == ddim) { 
                        // special case: P ~ 1/r => Integral ~ ln(r)
                        bond_len= pow(1.0*half_length+1.0,prng());
                    }
                    else { // normal case: P ~ 1/r^kappa with kappa!=1
                        bond_len= pow(lambda*prng()+1.0,dexponent);
                    }
                    // 2) Generate random d-dimensional versor
                    do{
                        bond_mod = 0.0;
                        for(dir=0;dir<(size_t)dims;dir++) {
                            vctr_coord[dir]=2*prng()-1.0;
                            bond_mod += pow(vctr_coord[dir], 2.0);
                        }
                    // Within unit sphere? - If not, discard!
                    } while (bond_mod > 1.0);
                    bond_mod=sqrt(bond_mod);
                    
                    // Build vector with random dir and length equal to bond_mod
                    for(dir=0;dir<dims;dir++){
                        vctr_coord[dir]=bond_len*vctr_coord[dir]/bond_mod;
                    } 
                    // Find site closest to this point. Recalculate length.
                    bond_mod=0.0;
                    ksite=0;
                    for(dir=0;dir<dims;dir++){
                        ksite += static_cast<size_t>(nearbyintl(vctr_coord[dir])*jump[dir]);
                        bond_mod+=nearbyintl(vctr_coord[dir])*nearbyintl(vctr_coord[dir]);
                    } 
                // make sure that bond is LONG
                } while (bond_len < pow(SHORT_BOND_LENGTH,2.0) ); 
            ksite=(isite+ksite+num_sites)%num_sites;
            break;
            } // End of long-bond method
    //	    break;
        }while(1); // repeat until a break is executed

	
	/* Final step: Modify neighbor list to include bond ij and count new bond */
	neighbor[isite][coordination[isite]]=ksite;
	neighbor[ksite][coordination[ksite]]=isite;
	coordination[isite]++;
	coordination[ksite]++;
	
    } /* END BOND LOOP */
    
    return;
}

// Code by Cristian Moukarzel, adopted by Chris Farrow
void SWNetwork::calculateBondWeights(void){

    size_t half_length = length/2;
    size_t isite,dir;
    size_t site_label,center_label,new_label;

    double bond_len;
    double total_sweight;

    // Calculate bond weights:
    // ----------------------
    //
    // Bond weights are defined as 1/(len)^{alpha}, where len is the length of
    // the bond, i.e. its distance to the origin. We are using helical BC's
    // here so there is a slight complication in order to define bond
    // lenghts. For a given site we will define its distance to the origin as
    // the minimum distance among all possible images of this site. This in
    // fact makes all distances smaller than approximately sqrt(d)*L/2.
    // 
    // In order to enforce this condition, the origin is defined to be at
    // {L/2,L/2,...,L/2} and bond lengths are calculated as distances to this
    // point (whose label is center_label).  

    // Calculate center_label
    center_label = 0;
    for(dir=0; dir < dims ;dir++)
    {
        center_label += half_length * jump[dir];
    }

    std::vector<size_t> coordinate(dims,0);

    // Calculate bond weights bond_weight[], and short-bond weight p_s
    total_sweight = 0.0;
    bond_weight[0] = 0.0;
    p_s = 0.0;
    for (isite=0; isite<num_sites ;isite++) {// loop over all possible links
        site_label = isite;
        // Label with respect to central point at {L/2,L/2,...,L/2}
        new_label = (isite - center_label + num_sites) % num_sites;
        bond_len = 0.0;
        // calculate d cartesian coordinates
        for(dir=0; dir < dims; dir++) {
            coordinate[dir] = site_label % length;
            site_label /= length;
            bond_len += (coordinate[dir]-half_length)*(coordinate[dir]-half_length);
        }

        if(new_label!=0) {// ignore central site
            bond_weight[new_label] = pow(bond_len,-alpha/2);
            // If bond is "short", add its weight to that of short bonds, and
            // store its identity.
            if(bond_len < pow(SHORT_BOND_LENGTH, 2.0) ){
                p_s += bond_weight[new_label];
                short_bond_label.push_back(new_label);
            }
            total_sweight+=bond_weight[new_label];
        }
    // end of loop over isites (i.e. over all possible bonds stemming from site
    // at center_label)
    }

    // Prepare the short_bond weight vector
    short_bond_weight.resize( short_bond_label.size() );

    // Normalize bond weights
    for (isite=0; isite<num_sites; isite++){
        bond_weight[isite] /= total_sweight;
    }
    p_s /= total_sweight;

    return;
}

size_t SWNetwork::getNumNeighbors(size_t site) {
    return coordination[site];
}

/* Note that the number of sites must be divisible by both the word_size
 * (defined in Lattice.h) and must be a power of length^dims.
 */
size_t scaleLength(size_t nsites, size_t dim) {
    size_t length = static_cast<size_t>(pow(nsites, 1.0/dim) );
    length -= length % word_size;
    nsites = static_cast<size_t>(pow(length, 1.0*dim));
    return nsites;
}
