#include <gtest/gtest.h>

#include "fk_mc.hpp"

#include "lattice/hypercubic.hpp"

using namespace fk;

int L = 8;
hypercubic_lattice<2> l1(L);

TEST(lattice, pos_to_index)
{

    for (int x=0; x<l1.msize(); x++) {
        EXPECT_EQ(x, l1.pos_to_index(l1.index_to_pos(x)));
        //std::cout << x << " -> " << l1.index_to_pos(x) << "->" << l1.pos_to_index(l1.index_to_pos(x)) << std::endl;
    };
}

TEST(lattice, hopping_m)
{
    fill_honeycomb(l1, 1.0);
    ///INFO(l1.hopping_m());
    // we only have 3 neighbors
    EXPECT_EQ(l1.hopping_m().sum(), -L*L*3);
}

