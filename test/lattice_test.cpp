#include <gtest/gtest.h>

#include "fk_mc.hpp"

#include "lattice/hypercubic.hpp"


using namespace fk;

size_t L = 8;
hypercubic_lattice<2> l1(L);

TEST(lattice, pos_to_index)
{

    for (int x=0; x<l1.msize(); x++) {
        EXPECT_EQ(x, l1.pos_to_index(l1.index_to_pos(x)));
        std::cout << x << " -> " << l1.index_to_pos(x) << "->" << l1.pos_to_index(l1.index_to_pos(x)) << std::endl;
    };
}

TEST(lattice, hopping_m)
{
    fill_nearest_neighbors(l1, 1.0);
    std::cout << l1.hopping_m() << std::endl;
}

TEST(lattice, bzpoints)
{
    hypercubic_lattice<2> t1(L);
    fill_triangular(t1, -1.0, -0.5);
    std::cout << l1.hopping_m() << std::endl;

    for (size_t i=0; i<t1.msize(); i++) {
        auto b = t1.get_bzpoint(i);
        FKDEBUG(b.ind_ << "->" << b << "<-" << t1.get_bzpoint(b.val_).ind_);
        ASSERT_EQ( t1.get_bzpoint(b.val_).ind_ , i);
    };

    auto bzpq = t1.get_bzpoint({{PI, PI}});
    FKDEBUG(bzpq.ind_ << "->" << bzpq << "<-" << t1.get_bzpoint(bzpq.val_).ind_);

    EXPECT_ANY_THROW({ auto bzpq1 = t1.get_bzpoint({{PI, PI+PI/7.}}); });

    auto bzpts = t1.get_all_bzpoints();
    for (auto x : bzpts) std::cout << x << " "; std::cout << std::endl;
}

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
