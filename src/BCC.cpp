/* BCC.h
 * 03/06/07 Created this as a subclass of LatticeRegular3d.
*/

#include "BCC.h"

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
void BCC<3>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( (i - x_val)/length ) % length;
    size_t z_val = (i - x_val - y_val*length) / length2;

    bool is_A = z_val % 2 == 0;
    int m = is_A ? -1 : 1;  // Multiply by -1 for A lattice
    short sx, sy;
    size_t x_new, y_new, z_new;

    for(int j=0; j<4; j++)
    {
        sx = m * B_shifts[j][0];
        sy = m * B_shifts[j][1];

        x_new = (x_val + sx + length) % length;
        y_new = (y_val + sy + length) % length;

        z_new = (z_val - 1 + length) % length;
        nbrs[2*j] = x_new + y_new * length + z_new * length2;

        z_new = (z_val + 1) % length;
        nbrs[2*j + 1] = x_new + y_new * length + z_new * length2;
    }
}

// PBC in X and Y directions. The lattice has a face.

template<>
size_t BCC<2>::getNumNeighbors(size_t i) {
    size_t x_val = i % length;
    size_t y_val = ( (i - x_val)/length ) % length;
    size_t z_val = (i - x_val - y_val*length) / length2;

    bool is_A = z_val % 2 == 0;
    bool botm_face = is_A && z_val == 0;
    bool top_face = !is_A && z_val == length - 1;

    if(botm_face || top_face) return 4;
    return 8;
}


template <>
void BCC<2>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( (i - x_val)/length ) % length;
    size_t z_val = (i - x_val - y_val*length) / length2;

    bool is_A = z_val % 2 == 0;
    int m = is_A ? -1 : 1;  // Multiply by -1 for A lattice
    short sx, sy;
    size_t x_new, y_new, z_new;
    
    bool botm_face = is_A && z_val == 0;
    bool top_face = !is_A && z_val == length - 1;

    nbrs.clear();
    for(int j=0; j<4; j++)
    {
        sx = m * B_shifts[j][0];
        sy = m * B_shifts[j][1];

        x_new = (x_val + sx + length) % length;
        y_new = (y_val + sy + length) % length;

        if(!botm_face)
        {
            z_new = (z_val - 1 + length) % length;
            nbrs.push_back(x_new + y_new * length + z_new * length2);
        }

        if(!top_face)
        {
            z_new = (z_val + 1) % length;
            nbrs.push_back(x_new + y_new * length + z_new * length2);
        }
    }
}


// PBC in X direction only. The structure has faces and edges.

template<>
size_t BCC<1>::getNumNeighbors(size_t i) {
    size_t x_val = i % length;
    size_t y_val = ( (i - x_val)/length ) % length;
    size_t z_val = (i - x_val - y_val*length) / length2;

    bool is_A = z_val % 2 == 0;
    bool botm_face = is_A && z_val == 0;
    bool top_face = !is_A && z_val == length - 1;
    bool front_face = is_A && y_val == 0;
    bool back_face = !is_A && y_val == length - 1;

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
void BCC<1>::setNbrs(size_t i) 
{
    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( (i - x_val)/length ) % length;
    size_t z_val = (i - x_val - y_val*length) / length2;

    bool is_A = z_val % 2 == 0;
    int m = is_A ? -1 : 1;  // Multiply by -1 for A lattice
    short sx, sy;
    size_t x_new, y_new, z_new;
    
    bool botm_face = is_A && z_val == 0;
    bool top_face = !is_A && z_val == length - 1;
    bool front_face = is_A && y_val == 0;
    bool back_face = !is_A && y_val == length - 1;

    nbrs.clear();

    for(int j=0; j<4; j++)
    {
        sx = m * B_shifts[j][0];
        sy = m * B_shifts[j][1];

        x_new = (x_val + sx + length) % length;
        y_new = (y_val + sy + length) % length;

        // Boundary conditions for A sublattice
        if(!botm_face && !(front_face && sy == -1))
        {
            z_new = (z_val - 1 + length) % length;
            nbrs.push_back(x_new + y_new * length + z_new * length2);
        }

        // Boundary conditions for B sublattice
        if(!top_face && !(back_face && sy == 1))
        {
            z_new = (z_val + 1) % length;
            nbrs.push_back(x_new + y_new * length + z_new * length2);
        }
    }
}


// No PBC. The lattice has faces, edges, and corners.

template<>
size_t BCC<0>::getNumNeighbors(size_t i) {
    size_t x_val = i % length;
    size_t y_val = ( (i - x_val)/length ) % length;
    size_t z_val = (i - x_val - y_val*length) / length2;

    bool is_A = z_val % 2 == 0;
    bool botm_face = is_A && z_val == 0;
    bool top_face = !is_A && z_val == length - 1;
    bool front_face = is_A && y_val == 0;
    bool back_face = !is_A && y_val == length - 1;
    bool left_face = is_A && x_val == 0;
    bool right_face = !is_A && x_val == length - 1;

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
void BCC<0>::setNbrs(size_t i) 
{

    /* Get the coordinates of the lattice site */
    size_t x_val = i % length;
    size_t y_val = ( (i - x_val)/length ) % length;
    size_t z_val = (i - x_val - y_val*length) / length2;

    bool is_A = z_val % 2 == 0;
    int m = is_A ? -1 : 1;  // Multiply by -1 for A lattice
    short sx, sy;
    size_t x_new, y_new, z_new;
    
    bool botm_face = is_A && z_val == 0;
    bool top_face = !is_A && z_val == length - 1;
    bool front_face = is_A && y_val == 0;
    bool back_face = !is_A && y_val == length - 1;
    bool left_face = is_A && x_val == 0;
    bool right_face = !is_A && x_val == length - 1;

    nbrs.clear();
    for(int j=0; j<4; j++)
    {
        sx = m * B_shifts[j][0];
        sy = m * B_shifts[j][1];

        x_new = (x_val + sx + length) % length;
        y_new = (y_val + sy + length) % length;

        // Boundary conditions for A sublattice
        if(!botm_face && !(front_face && sy == -1) && !(left_face && sx == -1))
        {
            z_new = (z_val - 1 + length) % length;
            nbrs.push_back(x_new + y_new * length + z_new * length2);
        }

        // Boundary conditions for B sublattice
        if(!top_face && !(back_face && sy == 1) && !(right_face && sx == 1))
        {
            z_new = (z_val + 1) % length;
            nbrs.push_back(x_new + y_new * length + z_new * length2);
        }
    }
}


// Let the compiler know we want these
template class BCC<0>;
template class BCC<1>;
template class BCC<2>;
template class BCC<3>;