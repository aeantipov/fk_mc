#include <gtest/gtest.h>

#include "fk_mc.hpp"

#include "lattice/hypercubic.hpp"
#include "lattice/chain.hpp"
#include "configuration.hpp"
#include "measures/polarization.hpp"


using namespace fk;


TEST(lattice, hopping)
{
    size_t L = 8;
    hypercubic_lattice<1> l1(L);
    chain_lattice l2(L);
    l1.fill(-0.8);
    l2.fill(-0.8, 0, 0);
    ASSERT_EQ(l2.msize(), L);
    ASSERT_EQ(l2.hopping_m().isApprox(l1.hopping_m()),true);
};

TEST(lattice, hopping2)
{
    size_t L = 8;
    chain_lattice l2(L);
    l2.fill(-0.8, 0.1, 0.07);
    ASSERT_EQ(l2.msize(), L);
    FKDEBUG(l2.hopping_m());
};

TEST(config, t1)
{
    size_t L = 36;
    chain_lattice l2(L);
    l2.fill(0.8, 0.0, 0.00);
    FKDEBUG(l2.hopping_m());

    double U = 2.0;
    double beta = 40;
    configuration_t config(l2, beta, U, -1.0 + U/2.*0, U/2.*0 );
    for (size_t x=0; x<L; x+=2)
       config.f_config_[x] = 0;
    config.calc_hamiltonian();
    FKDEBUG(config.hamilt_);
    config.calc_ed(true);
    FKDEBUG(config.ed_data().cached_spectrum.transpose());
    auto measure1 = measure_polarization<chain_lattice>(config,l2);
    measure1.accumulate(1.0);
};




int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
