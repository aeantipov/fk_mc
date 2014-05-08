#include "data_saveload.hpp"
#include <iostream>

namespace fk {

template <typename MC>
MC load_data(std::string output_file)
{
    H5::H5File input(output_file.c_str(),H5F_ACC_RDONLY);
    triqs::h5::group top(input);

    triqs::utility::parameters p;
    h5_read(top, "parameters", p);

    exit(0);
}

} // end of namespace fk
