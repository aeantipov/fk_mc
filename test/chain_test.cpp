#include <gtest/gtest.h>

#include "fk_mc.hpp"

#include "lattice/hypercubic.hpp"
#include "lattice/chain.hpp"


using namespace fk;


TEST(lattice, hopping)
{
    size_t L = 8;
    hypercubic_lattice<1> l1(L);
    chain_lattice l2(L);
    l1.fill(-0.8);
    l2.fill(-0.8, 0, 0);
    FKDEBUG(l1.hopping_m << " " << l2.hopping_m);
    ASSERT_EQ(l2.get_msize(), L);
    ASSERT_EQ(l2.hopping_m.isApprox(l1.hopping_m),true);
};

TEST(lattice, hopping2)
{
    size_t L = 8;
    chain_lattice l2(L);
    l2.fill(-0.8, 0.1, 0.07);
    ASSERT_EQ(l2.get_msize(), L);
    FKDEBUG(l2.hopping_m);
};




int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
