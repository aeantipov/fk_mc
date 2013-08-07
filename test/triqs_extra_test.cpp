#include <triqs/arrays.hpp>
#include "common.hpp"
#include "triqs_extra.hpp"

namespace tqa = triqs::arrays;


int main()
{
    std::vector<double> a(10,1.5);
    for (auto b: a) std::cout << b <<" " << std::flush;
    std::cout << std::endl;
    tqa::indexmaps::cuboid::domain_t<1> d1({int(a.size())});
    //tqa::indexmaps::cuboid::domain_t<1> d2(tqa::mini_vector<int, 1>({int(a.size())})); // FIXME!!
    tqa::indexmaps::cuboid::map<1,0,0> ind(d1);
    tqa::extra::weak_block<double> w(a.data(), a.size());
    tqa::extra::weak_view<double,1> e(ind, w);

    tqa::indexmaps::cuboid::domain_t<2> d2({2, 5});
    tqa::indexmaps::cuboid::map<2,0,0> ind2(d2);
    const tqa::extra::weak_view<double,2> e2(ind2, w);
    tqa::extra::weak_view<double,2> e3(ind2, w);

    std::cout << e << std::endl;
    std::cout << e2 << std::endl;

    e*=2.0;
    //e2*=2.0; // doesn't work as expected
    e3*=2.0;

    std::cout << e << std::endl;
    std::cout << e2 << std::endl;
    for (auto b: a) std::cout << b <<" " << std::flush;
    std::cout << std::endl;

    std::cout << fk::make_weak_view(a) << std::endl;
}
