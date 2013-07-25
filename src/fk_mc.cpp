#include "fk_mc.hpp"

namespace fk {

fk_mc::fk_mc(utility::parameters p)
{
   p.update(constructor_defaults()); 
}

} // end of namespace FK
