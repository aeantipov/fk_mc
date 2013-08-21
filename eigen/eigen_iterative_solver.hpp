#ifndef FK_EIGEN_ITERARIVE_INVERSE_SOLVER
#define FK_EIGEN_ITERARIVE_INVERSE_SOLVER

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Eigen{

// Iterative inverse eigensolver. See http://web.eecs.utk.edu/~dongarra/etemplates/node96.html
template<typename MatrixType, typename MatrixSolver=SimplicialLLT<MatrixType>>
class IterativeInverseEigenSolver {
public:

    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Index Index;
    typedef typename NumTraits<Scalar>::Real RealScalar;
    typedef typename internal::plain_col_type<MatrixType, RealScalar>::type RealVectorType;
    typedef Matrix<Scalar, Dynamic, Dynamic> DenseMatrixType;

    IterativeInverseEigenSolver()
     : m_eivec(),
       m_eivalues(),
       m_isInitialized(false),
       m_eigenvectorsOk(false),
       m_nbrConverged(0),
       m_nbrIterations(0)
    { }

    IterativeInverseEigenSolver(const MatrixType& A, Index nbrEigenvalues,
                                const RealVectorType &Eivalues_guess,
                                int options=ComputeEigenvectors, RealScalar tol=std::numeric_limits<RealScalar>::epsilon())
    {
        eigen_assert(nbrEigenvalues == Eivalues_guess.size());
        DenseMatrixType Evecs_guess(nbrEigenvalues, A.cols());
        Evecs_guess.setRandom(nbrEigenvalues, A.cols());
        compute(A, nbrEigenvalues, Eivalues_guess, Evecs_guess, options, tol);
    }

   IterativeInverseEigenSolver(const MatrixType& A, Index nbrEigenvalues,
                               const RealVectorType &Eivalues_guess, const DenseMatrixType &Evecs_guess,
                               int options=ComputeEigenvectors, RealScalar tol=std::numeric_limits<RealScalar>::epsilon())
   {
        eigen_assert(nbrEigenvalues == Eivalues_guess.size());
        eigen_assert(nbrEigenvalues == Evecs_guess.rows());
        eigen_assert(A.cols() == Evecs_guess.cols());
        compute(A, nbrEigenvalues, Eivalues_guess, Evecs_guess, options, tol);
   };

   void compute(const MatrixType& A, Index nbrEigenvalues, 
                const RealVectorType &Eivalues_guess, const DenseMatrixType& Evecs_guess, int options, 
                RealScalar tol = std::numeric_limits<RealScalar>::epsilon(), size_t max_n_iter = 1e8)
   {
        RealScalar diff = 1.0;
        for (size_t i=0; i<max_n_iter && diff > tol; ++i) {
        diff = 0.0;
        MY_DEBUG("!");
        }
   }; 

protected:
  Matrix<Scalar, Dynamic, Dynamic> m_eivec;
  Matrix<Scalar, Dynamic, 1> m_eivalues;
  ComputationInfo m_info;
  bool m_isInitialized;
  bool m_eigenvectorsOk;

  size_t m_nbrConverged;
  size_t m_nbrIterations;
};

} // end of namespace Eigen

#endif
