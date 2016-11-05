#include <boost/mpi/collectives.hpp>

#include "fk_mc/measures/eigenfunctions.hpp"

namespace fk {

measure_eigenfunctions::measure_eigenfunctions(configuration_t& in, std::vector<dense_m>& eigenfunctions):
    config(in), eigenfunctions_(eigenfunctions)
{ 
};

void measure_eigenfunctions::accumulate (double sign) 
{
    config.calc_ed(true);
    dense_m const& evecs = config.ed_data_.cached_evecs;
    eigenfunctions_.push_back(evecs);
    _Z++;
}

void measure_eigenfunctions::collect_results(boost::mpi::communicator const &c)
{
    int sum_Z;
    boost::mpi::reduce(c, _Z, sum_Z, std::plus<int>(), 0);
    std::vector<double> buffer;
    std::vector<dense_m> buffer2;
    size_t rows = eigenfunctions_[0].rows();
    size_t cols = eigenfunctions_[0].cols();
    size_t volume = rows*cols;
    for (size_t i=0; i<eigenfunctions_.size(); ++i) {
        buffer.resize(volume*c.size());
        boost::mpi::gather(c, eigenfunctions_[i].data(), volume, buffer, 0);
        for (int l=0; l<c.size(); ++l) { 
            buffer2.push_back(Eigen::Map<dense_m>(buffer.data() + l*volume, rows, cols));
            }
        };
    eigenfunctions_.swap(buffer2);
}


} // end of namespace FK
