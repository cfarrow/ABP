#ifndef NEIGHBORS_H
#define NEIGHBORS_H

/* Helper classes to simplify iteration over neighbors. 

Neighbor arrays used in LatticeRegular and LatticeRandom may over-allocate
memory for storing neighbor indices. Thus, the number of neighbors is explicitly
accounted for during iteration.

The Neighbors class can iterate over the passed neighbors array directly by
specifing safe=false. This should only be done if the neighbors array cannot
change while iterating.  When safe=true the neighbors array is copied before
iterating. If safety is too memory intensive, it is always possible to iterate
over neighbors using getNumNeihbors and getNbr.

*/

#include <algorithm>


class Neighbors {

    public:
        Neighbors(size_t* arr, const size_t& size, bool safe=true)
        : n(size), dealloc(safe) {
            if(safe){
                nbrs = new size_t [n];
                std::copy(arr, arr+n, nbrs);
            }
            else {
                nbrs = arr;
            }
        }

        ~Neighbors() {
            if(dealloc) delete [] nbrs;
        }


        auto begin() {return &nbrs[0];}
        auto end() {return begin() + n;}

    private:
        size_t* nbrs;
        bool dealloc;
        size_t n;
};


#endif