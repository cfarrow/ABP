/* BCC.h
 * 03/06/07 Created this as a subclass of LatticeRegular3d.
*/

#include "BCC.h"
#include "Point.h"

/* 
A BCC lattice comprises two intercalated cubic sub-lattices. We call them A and
B. Sub-lattice A has z values that are even and sub-lattice B has z values that
are odd. Any site on lattice A is at the center of a 8 sites on lattice B (it's
at the center of the cubic body, hence BCC).

If on sublattice 'A' (z even) the neighbors are:
x-1      y-1     z+/-1
x-1      y       z+/-1
x        y-1     z+/-1
x        y       z+/-1

If on sublattice 'B' (z odd) the neighbors are:
x        y       z+/-1
x        y+1     z+/-1
x+1      y       z+/-1
x+1      y+1     z+/-1

It's helpful later to define which sites are on a face for determining boundary
conditions.

x         | y         | z         | lattice
----------------------------------------------
0         | y         | even      | A
x         | 0         | even      | A
x         | y         | 0         | A
length-1  | y         | odd       | B
x         | length-1  | odd       | B
x         | y         | length-1  | B

A site can be on an edge or a corner if more than one condition is satisfied.
*/

// X and Y shifts on the B lattice to get to neighbors.
// The A shifts are the negative of these.
static const short B_shifts[4][2] = {
    {0, 0},
    {0, 1},
    {1, 0},
    {1, 1},
};

/* PBC in all directions */

template<>
size_t BCC<3>::getNumNeighbors(size_t) {
    return 8;
}


template <>
size_t BCC<3>::setNbrs(size_t i) 
{
    Point3d p{i, length};
    bool is_A = p.z % 2 == 0;
    int m = is_A ? -1 : 1;  // Multiply by -1 for A lattice
    short sx, sy;
    for(int j=0; j<4; j++)
    {
        sx = m * B_shifts[j][0];
        sy = m * B_shifts[j][1];
        nbrs[2*j]     = p.shift(sx, sy,  1);
        nbrs[2*j + 1] = p.shift(sx, sy, -1);
    }
    return 8;
}

// PBC in X and Y directions. The lattice has a face.

template<>
size_t BCC<2>::getNumNeighbors(size_t i) {
    Point3d p{i, length};
    bool is_A = p.z % 2 == 0;
    bool botm_face = is_A && p.z == 0;
    bool top_face = !is_A && p.z == b;

    if(botm_face || top_face) return 4;
    return 8;
}


template <>
size_t BCC<2>::setNbrs(size_t i) 
{
    Point3d p{i, length};
    size_t n{0};
    bool is_A = p.z % 2 == 0;
    int m = is_A ? -1 : 1;  // Multiply by -1 for A lattice
    short sx, sy;
    
    bool botm_face = is_A && p.z == 0;
    bool top_face = !is_A && p.z == b;

    for(int j=0; j<4; j++)
    {
        sx = m * B_shifts[j][0];
        sy = m * B_shifts[j][1];
        if(!botm_face)  nbrs[n++] = p.shift(sx, sy, -1);
        if(!top_face)   nbrs[n++] = p.shift(sx, sy,  1);
    }
    return n;
}


// PBC in X direction only. The structure has faces and edges.

template<>
size_t BCC<1>::getNumNeighbors(size_t i) {
    Point3d p{i, length};
    bool is_A = p.z % 2 == 0;
    bool botm_face = is_A && p.z == 0;
    bool top_face = !is_A && p.z == b;
    bool front_face = is_A && p.y == 0;
    bool back_face = !is_A && p.y == b;

    int num_faces = botm_face + front_face + top_face + back_face;

    switch (num_faces) {
        case 2:
            return 2;
            break;
        case 1:
            return 4;
            break;
        default:
            return 8;
            break;
    }
}


template<>
size_t BCC<1>::setNbrs(size_t i) 
{
    Point3d p{i, length};
    size_t n{0};
    bool is_A = p.z % 2 == 0;
    int m = is_A ? -1 : 1;  // Multiply by -1 for A lattice
    short sx, sy;
    
    bool botm_face = is_A && p.z == 0;
    bool top_face = !is_A && p.z == b;
    bool front_face = is_A && p.y == 0;
    bool back_face = !is_A && p.y == b;

    for(int j=0; j<4; j++)
    {
        sx = m * B_shifts[j][0];
        sy = m * B_shifts[j][1];
        // Boundary conditions for A sublattice
        if(!botm_face && !(front_face && sy == -1))
            nbrs[n++] = p.shift(sx, sy, -1);

        // Boundary conditions for B sublattice
        if(!top_face && !(back_face && sy == 1))
            nbrs[n++] = p.shift(sx, sy, 1);
    }
    return n;
}


// No PBC. The lattice has faces, edges, and corners.

template<>
size_t BCC<0>::getNumNeighbors(size_t i) {
    Point3d p{i, length};
    bool is_A = p.z % 2 == 0;
    bool botm_face = is_A && p.z == 0;
    bool top_face = !is_A && p.z == b;
    bool front_face = is_A && p.y == 0;
    bool back_face = !is_A && p.y == b;
    bool left_face = is_A && p.x == 0;
    bool right_face = !is_A && p.x == b;

    // A point can only be on a single sub-lattice, so we can sum the conditions
    // to determine if a point is on a corner, edge, face, or body.
    int num_faces = botm_face + front_face + left_face +
                    top_face + back_face + right_face;

    switch(num_faces){
        case 3:
            return 1;
            break;
        case 2:
            return 2;
            break;
        case 1:
            return 4;
            break;
        default:
            return 8;
            break;
    }
}


template<>
size_t BCC<0>::setNbrs(size_t i) 
{
    Point3d p{i, length};
    size_t n{0};
    bool is_A = p.z % 2 == 0;
    int m = is_A ? -1 : 1;  // Multiply by -1 for A lattice
    short sx, sy;
    
    bool botm_face = is_A && p.z == 0;
    bool top_face = !is_A && p.z == b;
    bool front_face = is_A && p.y == 0;
    bool back_face = !is_A && p.y == b;
    bool left_face = is_A && p.x == 0;
    bool right_face = !is_A && p.x == b;

    for(int j=0; j<4; j++)
    {
        sx = m * B_shifts[j][0];
        sy = m * B_shifts[j][1];

        // Boundary conditions for A sublattice
        if(!botm_face && !(front_face && sy == -1) && !(left_face && sx == -1))
            nbrs[n++] = p.shift(sx, sy, -1);

        // Boundary conditions for B sublattice
        if(!top_face && !(back_face && sy == 1) && !(right_face && sx == 1))
            nbrs[n++] = p.shift(sx, sy, 1);
    }
    return n;
}


// Let the compiler know we want these
template class BCC<0>;
template class BCC<1>;
template class BCC<2>;
template class BCC<3>;