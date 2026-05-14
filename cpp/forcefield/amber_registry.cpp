#include "amber_internal.hpp"

#include <mutex>
#include <shared_mutex>
#include <utility>

namespace xpongecpp {

std::unordered_map<std::string, ResidueType>& templates() {
    static std::unordered_map<std::string, ResidueType> registry;
    return registry;
}

std::unordered_map<std::string, Molecule>& molecule_templates() {
    static std::unordered_map<std::string, Molecule> registry;
    return registry;
}

std::shared_mutex& registry_mutex() {
    static std::shared_mutex mutex;
    return mutex;
}

void put_template(ResidueType residue_type) {
    configure_xponge_residue_links(residue_type);
    std::unique_lock lock(registry_mutex());
    molecule_templates().erase(residue_type.name());
    templates().insert_or_assign(residue_type.name(), std::move(residue_type));
}

}  // namespace xpongecpp
