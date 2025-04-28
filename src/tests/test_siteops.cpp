#include "../Lattice.h"

#include "gtest/gtest.h"

// There was an error in mask creation that was uncovered by switching to
// 64-bit from the 32-bit system on which the code was originally written and
// executed. Masks were creating by shifting a _integer_ literal '1' and then
// stored as a size_t. Shifts larger than the bit-width of an int resulted in
// no shift at all. (This has been fixed by explicityly creating a size_t{1} and
// then shifting.) As a reasult, half of all sites stored in the upper-half of a
// size_t word could not be modified by the code.  This test catches that error.

const size_t shift_3 = 0b1000;
const size_t shift_19 = 0b10000000000000000000;
const size_t shift_41 = 0b100000000000000000000000000000000000000000;


TEST(siteOpsTest, test_bit_mask)
{
    EXPECT_EQ(shift_3, bit_mask(3));

    if(word_size > 32) {
        EXPECT_EQ(shift_41, bit_mask(41));
    }
}


TEST(siteOpsTest, test_get_site_value)
{
    size_t* words = new size_t[8]{0, 0, 0, 0};
    size_t site_index;

    site_index = 2 * word_size + 3;
    EXPECT_EQ(false, get_site_value(site_index, words));

    words[2] |= shift_3;
    EXPECT_EQ(true, get_site_value(site_index, words));

    if(word_size > 32) {
        site_index = word_size + 41;
        EXPECT_EQ(false, get_site_value(site_index, words));

        words[1] |= shift_41;
        EXPECT_EQ(true, get_site_value(site_index, words));
    }

    delete [] words;
}

TEST(siteOpsTest, test_set_site_value)
{
    size_t* words = new size_t[8]{0, 0, 0, 0};
    size_t site_index;

    // Test flipping a single site.
    site_index = 2 * word_size + 3;
    set_site_value(site_index, true, words);
    EXPECT_EQ(0, words[0]);
    EXPECT_EQ(0, words[1]);
    EXPECT_EQ(shift_3, words[2]);
    EXPECT_EQ(0, words[3]);
    set_site_value(site_index, false, words);
    EXPECT_EQ(0, words[0]);
    EXPECT_EQ(0, words[1]);
    EXPECT_EQ(0, words[2]);
    EXPECT_EQ(0, words[3]);
 
    // Test flipping two sites
    set_site_value(site_index, true, words);
    set_site_value(2 * word_size + 19, true, words);
    EXPECT_EQ(0, words[0]);
    EXPECT_EQ(0, words[1]);
    EXPECT_EQ(shift_3 | shift_19, words[2]);
    EXPECT_EQ(0, words[3]);
    set_site_value(site_index, false, words);
    EXPECT_EQ(0, words[0]);
    EXPECT_EQ(0, words[1]);
    EXPECT_EQ(shift_19, words[2]);
    EXPECT_EQ(0, words[3]);

    if(word_size > 32) {
        // Test flipping a site beyond 32 bits
        site_index = word_size + 41;
        set_site_value(site_index, true, words);
        EXPECT_EQ(0, words[0]);
        EXPECT_EQ(shift_41, words[1]);
        EXPECT_EQ(shift_19, words[2]);
        EXPECT_EQ(0, words[3]);
        set_site_value(site_index, false, words);
        EXPECT_EQ(0, words[0]);
        EXPECT_EQ(0, words[1]);
        EXPECT_EQ(shift_19, words[2]);
        EXPECT_EQ(0, words[3]);
    }
}