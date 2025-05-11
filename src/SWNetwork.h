/* SWNetwork.h 11/10/07 - initial version of file
 *
 * Power law small-world network with poissonian distribution of bond lengths.
 *
 * Power-law networks are built by first taking a hypercubic lattice with no
 * bonds, and drawing bond-lengths l at random with probability P(l) ~
 * l^{-alpha). They are a generalization of SW networks with two control
 * parameters: the average coordination number gamma, and the power-law decay
 * exponent alpha. In the alpha->0 limit one recovers SW networks, and for large
 * alpha all bonds are short so the network is d-dimensional.
*/

#ifndef SWNET_H
#define SWNET_H

#include "Lattice.h"
#include "LatticeRandom.h"

size_t scaleLength(size_t nsites, size_t dim);

class SWNetwork: public LatticeRandom
{
    public:
        /*
         * nsites   --  The number of sites
         * id       --  An id number
         * dim      --  Number of dimensions to approach
         * alpha    --  Power law exponent
         * gamma    --  Average coordination number
         */
        SWNetwork(size_t nsites, size_t id, size_t dim = 2, double alpha = 2, double gamma = 6) : LatticeRandom(scaleLength(nsites, dim), id)
        {
            Setup(dim, alpha, gamma);
            generateBonds();
        }
        
        virtual ~SWNetwork() 
        {
            Cleanup();
        }
        
        virtual void generateBonds();
        virtual void activateSites();
        virtual size_t getNumNeighbors(size_t);
        size_t getDims() {return dim;}

    private:
        size_t dim;
        double alpha, gamma;

        size_t *coordination;
        double *bond_weight;
        std::vector<double> short_bond_weight;
        std::vector<size_t> short_bond_label;

        size_t *jump;

        void Setup(size_t dim, double alpha, double gamma);
        void Cleanup();
        
        void calculateBondWeights();

        // Constants
        size_t MAX_NUMBER_OF_SHORT_BONDS;
        double SHORT_BOND_LENGTH;
        double p_s; // Accumulated weight of all short bonds
};

#endif
