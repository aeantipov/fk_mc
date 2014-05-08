#pragma once
#include <triqs/h5.hpp>

#include "fk_mc.hpp"

namespace fk { 

// save data from solver to hdf5 file 
void save_data(const fk_mc& mc, triqs::utility::parameters p, std::string output_file, bool save_plaintext = false);

// construct a solver from given hdf5 file
fk_mc load_data(triqs::utility::parameters p, std::string output_file);

} // end of namespace fk

