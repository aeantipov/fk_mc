#include "fk_mc.hpp"

using namespace fk;

int main()
{
    square_lattice_traits<2> l1(4);

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

    l1.fill(-1.0);
    INFO(l1.get_hopping_matrix());

    triangular_lattice_traits t1(4);
    t1.fill(-1.0,-0.5);
    INFO(t1.get_hopping_matrix());

    return EXIT_SUCCESS;
}
