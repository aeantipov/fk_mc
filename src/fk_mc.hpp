#ifndef __FK_MC_HPP
#define __FK_MC_HPP

#include "common.hpp"
//#include "parameters.hpp"
#include <boost/mpi/communicator.hpp>

using namespace triqs;

namespace fk {

 class fk_mc {

  // define the communicator, here MPI_COMM_WORLD
  boost::mpi::communicator world;

  utility::parameter_defaults constructor_defaults() const;
  utility::parameter_defaults solve_defaults() const;

  public:

  fk_mc(utility::parameters p);
  void solve(utility::parameters p);

 };

}; // end of namespace FK

#endif // #ifndef __FK_MC_HPP
