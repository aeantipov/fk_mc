#include <gtest/gtest.h>

#include "fk_mc.hpp"

#include "lattice/hypercubic.hpp"
#include "measures/stiffness.hpp"

inline void print_section (const std::string& str) { std::cout << std::string(str.size(),'-') << "\n" << str << std::endl; }

using namespace fk;

void test_stiffness(double U, double beta, std::vector<int> f_config, double compare_val)
{ 
    boost::mpi::communicator comm;
    typedef hypercubic_lattice<2> lattice_t;
    int L = std::lround(sqrt(f_config.size()));
    std::cout << "L = " << L << std::endl;
    double mu = U/2;
    double e_f = 0.0;
    double t = 1.0;

    lattice_t lattice(L);
    lattice.fill(t);
    configuration_t config(lattice, beta, U, mu, mu+e_f);
    Eigen::ArrayXi my_config(lattice.get_msize()); 
    std::copy(f_config.begin(), f_config.end(), &my_config[0]);

    config.f_config_ = my_config;
    std::cout << "f-electron config : " << config.f_config_.transpose() << std::endl;

    config.calc_hamiltonian();
    config.calc_ed(true);

    std::cout << "Spectrum : " << config.ed_data().cached_spectrum.transpose() << std::endl;
    
    std::vector<double> stiffness_vals;
    std::vector<double> wgrid = {0.0};
    std::vector<std::vector<double>> cond_history;
    auto st_m = measure_stiffness<lattice_t>(config,lattice,stiffness_vals,cond_history,wgrid,0.05);

    st_m.accumulate(1.0);
    st_m.collect_results(comm);
    for (auto& v:stiffness_vals) std::cout << v << " " << std::flush; std::cout << std::endl;

    double obs, obs_err, comp;
    std::string title;
    #define compare_me(check) {                                                 \
        print_section(title);                                                   \
        std::cout << obs << " +/- " << obs_err  << " == " << comp << std::endl; \
        if (check) EXPECT_NEAR(obs, comp, obs_err);                             \
        std::cout << std::endl; };

    title = "stiffness";
    obs = stiffness_vals[0];
    obs_err = 1e-3;
    comp = compare_val; 
    compare_me(true);
}

TEST(stiffness, test00) { test_stiffness(0.0, 1000.0, {0,1,0,1,1,0,1,0,0,1,0,0,1,0,1,0,1,1,1,0,1,1,0,1,1}, 1.31597); };
TEST(stiffness, test01) { test_stiffness(2.0, 1000.0, {0,1,0,1,1,0,1,0,0,1,0,0,1,0,1,0,1,1,1,0,1,1,0,1,1}, 0.287547); };
TEST(stiffness, test02) { test_stiffness(6.0, 1000.0, {0,1,0,1,1,0,1,0,0,1,0,0,1,0,1,0,1,1,1,0,1,1,0,1,1}, 0.00260831); };
TEST(stiffness, test10) { test_stiffness(0.37, 1000.0, {0,1,0,0,0,0,1,1,0,1,1,0,0,1,1,0,1,1,0,1,1,0,1,0,0,0,1,1,0,1,1,1,1,1,0,1,1,0,1,1,1,1,0,0,0,1,0,1,1}, 1.26046); }

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
