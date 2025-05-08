#ifndef NEIGHBORS_H
#define NEIGHBORS_H


/* Helper class to simplify iteration over neighbors. 

Neighbor arrays used in LatticeRegular and LatticeRegular may over-allocate
memory for storing neighbor indices. Thus, the number of neighbors is explicitly
accounted for during iteration.

*/
class Neighbors {

    public:
        Neighbors(const size_t* arr, const size_t& size) : nbrs(arr), n(size) {}

        auto begin() {return &nbrs[0];}
        auto end() {return begin() + n;}

    private:
        const size_t* nbrs;
        const size_t n;
};

#endif