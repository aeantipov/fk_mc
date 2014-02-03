#include "fk_mc.hpp"

using namespace fk;

int main()
{
    hypercubic_lattice<2> l1(4);

    MY_DEBUG(l1.pos_to_index({0,0}));
    MY_DEBUG(l1.pos_to_index({0,2}));
    MY_DEBUG(l1.pos_to_index({1,0}));
    MY_DEBUG(l1.pos_to_index({2,1}));

    MY_DEBUG(l1.index_to_pos(0));
    MY_DEBUG(l1.index_to_pos(2));
    MY_DEBUG(l1.index_to_pos(4));
    MY_DEBUG(l1.index_to_pos(9));

    bool success;
    std::array<size_t, 2> pos1 = { 2, 3};
    success = (l1.index_to_pos(l1.pos_to_index(pos1)) == pos1);
    INFO(pos1 << "==" << l1.index_to_pos(l1.pos_to_index(pos1)) << " = " << std::boolalpha << success);
    if (!success) return EXIT_FAILURE;
    success = (l1.pos_to_index(l1.index_to_pos(11)) == 11);
    INFO(11 << "==" << l1.pos_to_index(l1.index_to_pos(11)) << " = " << std::boolalpha << success);
    if (!success) return EXIT_FAILURE;

    try {
    l1.fill(-1.0);
        }
    catch (std::exception &e){MY_DEBUG(e.what());};
    INFO(l1.hopping_m);

    Eigen::VectorXcd a1(l1.get_msize()); 
    a1.setZero();
    a1[3]=1.;
    a1[8]=1.;
    MY_DEBUG(a1.transpose());
    auto af1 = l1.FFT(a1, FFTW_FORWARD);
    MY_DEBUG(af1.transpose());
    auto a2 = l1.FFT(af1, FFTW_BACKWARD);
    MY_DEBUG(a2.transpose())
    if (a2 != a1) return EXIT_FAILURE;

    triangular_lattice t1(4);
    t1.fill(-1.0,-0.5);
    INFO(t1.hopping_m);

    return EXIT_SUCCESS;
}
