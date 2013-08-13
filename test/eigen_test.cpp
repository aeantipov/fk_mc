#include "common.hpp"
#include <chrono>
#include <unsupported/Eigen/ArpackSupport>
#include <triqs/mc_tools/random_generator.hpp>


using namespace fk;

int main(int argc, char* argv[])
{

    typedef Eigen::SparseMatrix<double> sparse_m;
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dense_m;
    triqs::mc_tools::random_generator RNG("mt19937", 23432);

    size_t size = 10;
    size_t nonzero_elems = size*3;
    //EMatrixType<double> a(10,10);
    //a.setRandom();
    sparse_m a(size,size); 
    a.reserve(2*nonzero_elems);
    for (size_t l=0; l<nonzero_elems; l++) { 
        size_t i = RNG(size), j = i+RNG(size-i); 
        while(a.coeffRef(i,j)!=0) {i = RNG(size); j = RNG(size);};
        a.coeffRef(i,j) = RNG(2.0)-1.;
        a.coeffRef(j,i) = a.coeffRef(i,j);
        };
    INFO(a);

    dense_m b(a);
    INFO(b);
    Eigen::SelfAdjointEigenSolver<dense_m> s(b);
    INFO(s.eigenvalues());
    //a = a.selfadjointView<Eigen::Upper>();
    INFO("");
    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<sparse_m> solver(a,2,"SA",Eigen::EigenvaluesOnly,1e-5);
    INFO(solver.eigenvalues());
    return EXIT_SUCCESS;
}
