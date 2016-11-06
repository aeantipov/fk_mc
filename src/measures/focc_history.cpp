#include <boost/mpi/collectives.hpp>

#include "fk_mc/measures/focc_history.hpp"

namespace fk {

void measure_focc::measure (double sign) 
{
    const auto& focc_current = config.f_config_;
    for (size_t i=0; i<focc_current.size(); ++i) { focc_[i].push_back(focc_current(i)); };
    _Z++;
}

void measure_focc::collect_results(boost::mpi::communicator const &c)
{
    int sum_Z;
    boost::mpi::reduce(c, _Z, sum_Z, std::plus<int>(), 0);
    std::vector<double> focc_tmp;
    for (size_t i=0; i<focc_.size(); ++i) {
        focc_tmp.resize(focc_[i].size()*c.size());
        boost::mpi::gather(c, focc_[i].data(), focc_[i].size(), focc_tmp, 0);
        focc_[i].swap(focc_tmp);
        };
}

} // end of namespace FK
