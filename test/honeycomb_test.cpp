#include <gtest/gtest.h>

#include "fk_mc.hpp"

#include "lattice/honeycomb.hpp"

using namespace fk;

int L = 8;
honeycomb_lattice l1(L);

TEST(lattice, pos_to_index)
{

    for (int x=0; x<l1.msize(); x++) {
        EXPECT_EQ(x, l1.pos_to_index(l1.index_to_pos(x)));
        //std::cout << x << " -> " << l1.index_to_pos(x) << "->" << l1.pos_to_index(l1.index_to_pos(x)) << std::endl;
    };
}

TEST(lattice, hopping_m)
{
    l1.fill(1.0);
    ///INFO(l1.hopping_m());
    // we only have 3 neighbors
    EXPECT_EQ(l1.hopping_m().sum(), -L*L*3);
}
/*
TEST(lattice, FFT)
{
    honeycomb_lattice t1(L);
    t1.fill(1.0);

    Eigen::MatrixXcd hop_m = t1.hopping_m().cast<std::complex<double>>(); 
    Eigen::MatrixXcd disp = t1.FFT(hop_m.row(0), FFTW_FORWARD);

    std::cout << disp.real().cast<double>() << std::endl;
}
*/
/*
TEST(lattice, fft)
{
    Eigen::ArrayXcd a1(l1.msize());
    a1.setZero();
    a1[3]=1.;
    a1[8]=1.;
    FKDEBUG(a1.transpose());
    auto af1 = l1.FFT(a1, FFTW_FORWARD);
    FKDEBUG(af1.transpose());
    auto a2 = l1.FFT(af1, FFTW_BACKWARD);
    FKDEBUG(a2.transpose())
    EXPECT_EQ(a2.isApprox(a1), true);
}


TEST(lattice, fft2)
{ 
    Eigen::ArrayXi a1(l1.msize());
    a1.setZero();
    for (int i=0; i<l1.msize(); i++) {
        auto pos = l1.index_to_pos(i);
        int s = -((pos[0] + pos[1])%2 * 2 - 1);
        a1[i] = s;
        };
    std::cout << a1.transpose() << " == " << std::endl << l1.ft_pi_array_.transpose() << std::endl;
    ASSERT_EQ(l1.ft_pi_array_.isApprox(a1),true);
};



int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
*/
