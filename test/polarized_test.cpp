#include <gtest/gtest.h>

#include "fk_mc.hpp"

#include "lattice/polarized.hpp"


using namespace fk;

size_t L = 4;
polarized_lattice<2> l1(L);

    double U = 4.0;
    double mu = U/2;
    double e_f = 0.0;
    double t = 1.0;
    double beta = 5.0;



TEST(lattice, pos_to_index)
{

    for (int x=0; x<l1.get_msize(); x++) {
        EXPECT_EQ(x, l1.pos_to_index(l1.index_to_pos(x)));
        std::cout << x << " -> " << l1.index_to_pos(x) << "->" << l1.pos_to_index(l1.index_to_pos(x)) << std::endl;
    };
}

TEST(lattice, hopping_m)
{
    l1.fill(1.0, 10.0);
    INFO(l1.hopping_m);
}

TEST(lattice, spectrum)
{
 //   l1.fill(1.0, 10.0);
    configuration_t config(l1, beta, U, mu, mu+e_f);
    // checkerboard
    Eigen::ArrayXi my_config(l1.get_msize()); my_config.setZero();
    for (size_t x=0; x<L; x+=1)
        for (size_t y=0; y<L; y+=2) {
            my_config(l1.pos_to_index({x,y+x%2}))=1;
        }
    config.f_config_ = my_config;
    std::cout << "f-electron config : " << config.f_config_.transpose() << std::endl;

    config.calc_hamiltonian();
    config.calc_ed(true);

    std::cout << "Spectrum : " << config.ed_data().cached_spectrum.transpose() << std::endl;
    FKDEBUG(config.ed_data().cached_evecs);

}


int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
