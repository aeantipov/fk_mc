#include <gtest/gtest.h>
#include "lattice/graph.hpp"

using namespace fk;


TEST(lattice, construct)
{
    int ndim = 1;
    typedef std::vector<int> pos_t;
    std::vector<pos_t> sites;


    int Nx = 4, Ny = 4;
    int npts = Nx * Ny;
    int p = 0;
    for (int x = 0; x < Nx; ++x) {
        for (int y = 0; y < Ny; ++y) {
            sites.push_back(pos_t({x, y}));
        }
    }

    for (auto p : sites) { std::cout << "(" << p[0] << ";" << p[1] << ")" << std::flush;}
    std::cout << std::endl;
    typename abstract_lattice::sparse_m hop_m(Nx*Ny, Nx*Ny);
    //hop_m.coeffRef(p1, p2) = m;

    graph_lattice l(sites);

    std::cout << "Number of dimensions = " << l.ndim() << std::endl;
    ASSERT_EQ(l.ndim(), 2);
}

TEST(lattice, read_hdf5)
{
    alps::hdf5::archive ar("graph_lattice_test.h5", "r");
    read_lattice_hdf5(ar, "lattice");
    ar.close();
}

