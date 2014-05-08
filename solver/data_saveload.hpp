#pragma once
#include <triqs/h5.hpp>

#include "fk_mc.hpp"

namespace fk { 

// save data from solver to hdf5 file 
template <typename MC>
void save_data(const MC& mc, triqs::utility::parameters p, std::string output_file, bool save_plaintext = false);

// construct a solver from given hdf5 file. The observables are populated with existing data
template <typename MC>
MC load_data(std::string output_file);

} // end of namespace fk

#include "data_save.hxx"

