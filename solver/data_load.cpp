#include <fstream>

#include "data_saveload.hpp"

namespace fk {

fk_mc load_data(triqs::utility::parameters p, std::string output_file);

void print_section (const std::string& str)
{
  std::cout << std::string(str.size(),'=') << std::endl;
  std::cout << str << std::endl;
}

} // end of namespace fk
