#include "bindings_internal.hpp"

namespace xpongecpp {

PYBIND11_MODULE(_core, m) {
    bind_core_module(m);
    bind_io_module(m);
    bind_forcefield_module(m);
    bind_assign_module(m);
}

}  // namespace xpongecpp
