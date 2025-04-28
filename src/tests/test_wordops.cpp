#include "../Lattice.h"

#include "gtest/gtest.h"

// There was an error in mask creation that was uncovered by switching to
// 64-bit from the 32-bit system on which the code was originally written and
// executed. Masks were creating by shifting a _integer_ literal '1' and then
// stored as a size_t. Shifts larger than the bit-width of an int resulted in
// no shift at all. (This has been fixed by explicityly creating a size_t{1} and
// then shifting.) As a reasult, half of all sites stored in the upper-half of a
// size_t word could not be modified by the code.  This test catches that error.
TEST(wordopsTest, testbitmask)
{
    size_t shift, answer;

    shift = 1;
    answer = 0b10;
    ASSERT_EQ(answer, bit_mask(shift));

    if(sizeof(size_t) > 4) {
        shift = 41;
        answer = 0b100000000000000000000000000000000000000000;
        ASSERT_EQ(answer, bit_mask(shift));
    }
}