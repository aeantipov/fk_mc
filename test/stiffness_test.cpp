#include <gtest/gtest.h>

#include "fk_mc.hpp"

#include "lattice/hypercubic.hpp"
#include "measures/stiffness.hpp"


using namespace fk;

size_t L = 4;
typedef hypercubic_lattice<2> lattice_t;
lattice_t lattice(L);

double U = 1.0;
double mu = U/2;
double e_f = 0.0;
double t = 1.0;
double beta = 100.0;

TEST(lattice, spectrum)
{
    boost::mpi::communicator comm;
 //   lattice.fill(1.0, 10.0);
    lattice.fill(1.0);
    configuration_t config(lattice, beta, U, mu, mu+e_f);
    // checkerboard
    Eigen::ArrayXi my_config(lattice.get_msize()); my_config.setZero();

    Eigen::Matrix4i positions;
    positions << 
         1, 0, 1, 0,
         1, 1, 0, 0,
         1, 0, 0, 1,
         1, 0, 1, 0
        ;
    for (int x=0; x<L; x++)
        for (int y=0; y<L; y++) {
            my_config(lattice.pos_to_index({{x,y}}))=positions(x,y);
        }
    config.f_config_ = my_config;
    std::cout << "f-electron config : " << config.f_config_.transpose() << std::endl;

    config.calc_hamiltonian();
    config.calc_ed(true);

    std::cout << "Spectrum : " << config.ed_data().cached_spectrum.transpose() << std::endl;
    FKDEBUG(config.ed_data().cached_evecs);
    
    std::vector<double> stiffness_vals;
    auto st_m = measure_stiffness<lattice_t>(config,lattice,stiffness_vals);

    st_m.accumulate(1.0);
    st_m.collect_results(comm);
    for (auto& v:stiffness_vals) std::cout << v << " " << std::flush; std::cout << std::endl;
}

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
