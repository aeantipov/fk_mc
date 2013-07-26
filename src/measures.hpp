#ifndef __FK_MC_MEASURES_HPP_
#define __FK_MC_MEASURES_HPP_

#include "common.hpp"


namespace fk {

// a measurement: the magnetization
struct dummy_measure {

  dummy_measure() = default;
 
  void accumulate(double sign) { }

  void collect_results(boost::mpi::communicator const &c) {

    double sum_Z, sum_M;
//    boost::mpi::reduce(c, Z, sum_Z, std::plus<double>(), 0);
//    boost::mpi::reduce(c, M, sum_M, std::plus<double>(), 0);

    if (c.rank() == 0) {
//      std::cout << "Magnetization: " << sum_M / sum_Z << std::endl << std::endl;
    }

  }

};

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURES_HPP_
