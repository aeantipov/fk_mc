#include <gtest/gtest.h>

#include "lattice/hypercubic.hpp"
#include "fk_mc.hpp"
#include <boost/mpi/environment.hpp>
#include <chrono>

using namespace fk;

TEST(config, test1)
{
    size_t L = 4;

    double U = 1.0;
    double mu = U/2;
    double e_f = 0.0;
    double t = -1.0;
    double T = 0.15;
    double beta = 1.0/T;

    hypercubic_lattice<2> lattice(L);
    fill_nearest_neighbors(lattice, -1.0);

    FKDEBUG(lattice.hopping_m());

    configuration_t config(lattice, 1.0, U, mu, mu+e_f);
    configuration_t config2(config);
    configuration_t config3(configuration_t(lattice, 1.0, U, mu, mu+e_f));
    config3 = config2;
    config3 = configuration_t(lattice, 1.0, U, mu, mu+e_f);

    Eigen::ArrayXi my_config(lattice.msize()); my_config.setZero();
    for (int x=0; x<L; x+=1)
        for (int y=0; y<L; y+=2) {
            my_config(lattice.pos_to_index({{x,y+x%2}}))=1;
        }
    config.f_config_ = my_config;
    FKDEBUG(config.f_config_.transpose());
    config.calc_hamiltonian();
    config.calc_ed();
    FKDEBUG(config.ed_data().cached_spectrum.transpose());
    //FKDEBUG(config.cached_evecs);
}

// test ff-energy
TEST(config, ff_energy)
{
    size_t L = 12;
    double U = 0.0;
    double mu = U/2;
    double e_f = 0.0;

    std::vector<double> W = { 0, 1, 2};

    hypercubic_lattice<1> lattice(L);
    fill_nearest_neighbors(lattice, -1.0);
    configuration_t config(lattice, 1.0, U, mu, mu+e_f, W);

    Eigen::ArrayXi my_config(lattice.msize()); my_config.setZero();
    for (int x=0; x<L; x+=1) {
            my_config(x) = x % 2;
        }
    my_config(0) = 1;
    config.f_config_ = my_config;

    double e_ff = config.calc_ff_energy(); 
    double e_ff_comp = 2 * L/2 * W[2] + 4 * W[1];

    std::cout << "configuration = " << my_config.transpose() << std::endl;
    std::cout << "ff energy = " << e_ff << " == " << e_ff_comp << std::endl;
    ASSERT_NEAR(e_ff, e_ff_comp, 1e-15);
}

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
