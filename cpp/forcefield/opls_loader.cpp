#include "nonamber_internal.hpp"

namespace xpongecpp {

Molecule load_opls_itp_file(const std::filesystem::path& filename) {
    return load_gromacs_topology_file(filename);
}

}  // namespace xpongecpp
