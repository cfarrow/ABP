#ifndef ABPPROCS_H
#define ABPPROCS_H

#include <queue>

#include "Lattice.h"
#include "rand.h"

/* Class for holding graph info. */
class GraphInfo
{
    public:
        const static size_t  HEX_LATTICE = 0; // hexagonal
        const static size_t  SQU_LATTICE = 1; // square
        const static size_t  TRI_LATTICE = 2; // triangular
        const static size_t  JAK_LATTICE = 3; // union jack (2-d, 8-coordinated)
        const static size_t  CUB_LATTICE = 4; // cubic
        const static size_t  BCC_LATTICE = 5; // bcc
        const static size_t  C4D_LATTICE = 6; // 4D cubic
        const static size_t  FIZ_LATTICE = 7; // Fixed coordination random
        const static size_t  SQR_LATTICE = 8; // Square-Random
        const static size_t  TRR_LATTICE = 9; // Triangle-Random
        const static size_t  CUR_LATTICE = 10; // Cubic-Random
        const static size_t  SMW_LATTICE = 11; // Small-world

        GraphInfo() { 
            type = coord = dim = 0;
            pbc = 0;
            fb = alpha = gamma = 0;
        }

        size_t type; // Lattice type
        size_t coord; // coordination
        short pbc; // periodic boundaries
        double fb; // fraction of long-range bonds in reg/rand graph
        size_t dim; // Dimension of small world network
        double alpha; // SW exponent
        double gamma; // SW average coordination

        // Methods!
        void promptGraphType();
        void printType( std::ostream &ostr);
        void printGraphInfo( std::ostream &ostr, size_t &nsites, size_t &mcull, size_t &num_samples);
        Lattice* getLattice(size_t len);
};


size_t removeActiveSite( Lattice *L, std::queue< size_t > &Q, size_t &t, size_t &k, bounded_rng_type& rng);
void collectToQueue( Lattice *L, std::queue< size_t > &Q, size_t cull );
size_t cullSites( Lattice *L, std::queue< size_t > &Q, size_t cull );
void promptLength( size_t &len );
void promptCullParam( size_t &cull);
void promptNumSamples( size_t &num_samples);
void processCommandLine(int argc, char **argv, size_t &len, size_t &mcull, size_t &num_samples, GraphInfo &gi);

#endif
