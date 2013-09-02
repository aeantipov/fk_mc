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
                                int options=ComputeEigenvectors, RealScalar tol=std::numeric_limits<RealScalar>::epsilon(),
                                size_t max_n_iter = 1e8
                                )
    {
        eigen_assert(nbrEigenvalues == Eivalues_guess.size());
        DenseMatrixType Evecs_guess;
        Evecs_guess.setRandom(A.rows(), nbrEigenvalues);
        compute(A, nbrEigenvalues, Eivalues_guess, Evecs_guess, options, tol, max_n_iter);
    }

   IterativeInverseEigenSolver(const MatrixType& A, Index nbrEigenvalues,
                               const RealVectorType &Eivalues_guess, const DenseMatrixType &Evecs_guess,
                               int options=ComputeEigenvectors, RealScalar tol=std::numeric_limits<RealScalar>::epsilon(), 
                               size_t max_n_iter = 1e8
                              )
   {
        eigen_assert(nbrEigenvalues == Eivalues_guess.size());
        eigen_assert(nbrEigenvalues == Evecs_guess.cols());
        eigen_assert(A.rows() == Evecs_guess.rows());
        compute(A, nbrEigenvalues, Eivalues_guess, Evecs_guess, options, tol, max_n_iter);
   };

   void compute(const MatrixType& A, Index nbrEigenvalues, 
                const RealVectorType &Eivalues_guess, const DenseMatrixType& Evecs_guess, int options, 
                RealScalar tol = std::numeric_limits<RealScalar>::epsilon(), size_t max_n_iter = 1e8)
   {
        MatrixSolver OP;
        // See algorithm at http://web.eecs.utk.edu/~dongarra/etemplates/node96.html
        size_t n_evals = Eivalues_guess.size();
        for (size_t n=0; n < n_evals; ++n) { 
            RealScalar sigma = Eivalues_guess(n);

            RealScalar theta;
            MatrixType A2(A);
            for (size_t j=0; j<A2.rows(); j++) A2.coeffRef(j,j) -= sigma;
            OP.compute(A2);
            //MY_DEBUG(OP.determinant());
            RealVectorType y = Evecs_guess.col(n); y/=y.norm();
            RealVectorType v(y);

            RealScalar diff = 1.0;
            for (size_t i=0; i<max_n_iter && diff > tol; ++i) {
                v = y/y.norm();
                y = A2.ldlt().solve(v); //OP.solve(v);
                theta = v.transpose()*y;
                diff = (y - theta*v).norm();
                MY_DEBUG("theta = " << theta << ";diff = " << diff);
                }
            MY_DEBUG("sigma = " << sigma << " theta = " << theta);
            MY_DEBUG("Eval: " << sigma + 1.0/theta);
            MY_DEBUG("Input guess: " << Evecs_guess.col(n).transpose());
            MY_DEBUG("Output evec: " << y.transpose()/theta);
            MY_DEBUG((A*Evecs_guess.col(n)/sigma).transpose());
            MY_DEBUG((A*Evecs_guess.col(n)/(sigma + 1.0/theta)).transpose());
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
