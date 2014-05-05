#include "fk_mc.hpp"

using namespace fk;

int main()
{
    size_t L = 4;
    hypercubic_lattice<2> l1(L);

    DEBUG(l1.pos_to_index({0,0}));
    DEBUG(l1.pos_to_index({0,2}));
    DEBUG(l1.pos_to_index({1,0}));
    DEBUG(l1.pos_to_index({2,1}));

    DEBUG(l1.index_to_pos(0));
    DEBUG(l1.index_to_pos(2));
    DEBUG(l1.index_to_pos(4));
    DEBUG(l1.index_to_pos(9));

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
    catch (std::exception &e){DEBUG(e.what());};
    INFO(l1.hopping_m);

    Eigen::ArrayXcd a1(l1.get_msize()); 
    a1.setZero();
    a1[3]=1.;
    a1[8]=1.;
    DEBUG(a1.transpose());
    auto af1 = l1.FFT(a1, FFTW_FORWARD);
    DEBUG(af1.transpose());
    auto a2 = l1.FFT(af1, FFTW_BACKWARD);
    DEBUG(a2.transpose())
    if (!a2.isApprox(a1)) return EXIT_FAILURE;

    triangular_lattice t1(L);
    t1.fill(-1.0,-0.5);
    INFO(t1.hopping_m);


    for (size_t i=0; i<L*L; i++) {
        auto b = t1.get_bzpoint(i);
        DEBUG(b.ind_ << "->" << b << "<-" << t1.get_bzpoint(b.val_).ind_);
        if ( t1.get_bzpoint(b.val_).ind_ != i) return EXIT_FAILURE;
    };

    auto bzpq = t1.get_bzpoint({PI, PI});
    DEBUG(bzpq.ind_ << "->" << bzpq << "<-" << t1.get_bzpoint(bzpq.val_).ind_);

    try { auto bzpq1 = t1.get_bzpoint({PI, PI+PI/7.}); }
    catch(triqs::runtime_error const & e) { std::cout  << "Caught exception "<< e.what() << std::endl;}

    auto bzpts = t1.get_all_bzpoints();
    for (auto x : bzpts) std::cout << x << " "; std::cout << std::endl;

    return EXIT_SUCCESS;
}
