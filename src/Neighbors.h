#ifndef NEIGHBORS_H
#define NEIGHBORS_H


/* Helper class to simplify iteration over neighbors. 

Neighbor arrays used in LatticeRegular and LatticeRegular may over-allocate
memory for storing neighbor indices. Thus, the number of neighbors is explicitly
accounted for during iteration.

*/
class Neighbors {

    public:
        Neighbors(const size_t* arr, const size_t& size) : n(size) {
            // We have to make a copy because other functions can change the
            // array while we're iterating.
            nbrs = new size_t [n];
            std::copy(arr, arr + n, nbrs);
        }

        ~Neighbors(){
            delete [] nbrs;
        }

        auto begin() {return &nbrs[0];}
        auto end() {return begin() + n;}

    private:
        size_t* nbrs;
        const size_t n;
};

#endif