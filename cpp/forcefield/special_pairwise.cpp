#include "nonamber_internal.hpp"

#include <fstream>
#include <optional>
#include <stdexcept>

namespace xpongecpp {

void load_sw_parameter_file(const std::filesystem::path& filename, Molecule& molecule) {
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open SW parameter file: " + filename.string());
    }
    for (auto& atom : molecule.atoms) {
        atom.sw_type = atom.type;
        if (atom_type_registry().find(atom.type) == atom_type_registry().end()) {
            register_external_lj(atom.type, 0.0, 0.0);
        }
    }
    for (std::string line; std::getline(input, line);) {
        const auto words = split_ws(strip_comment(line));
        if (words.size() < 11) {
            continue;
        }
        const auto name = words[0];
        molecule.add_sw_type(name, std::stod(words[7]), std::stod(words[8]), std::stod(words[1]),
                             std::stod(words[9]), std::stod(words[10]), std::stod(words[3]),
                             std::stod(words[5]), std::stod(words[2]), std::stod(words[4]), std::stod(words[6]));
    }
}

void load_edip_parameter_file(const std::filesystem::path& filename, Molecule& molecule) {
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open EDIP parameter file: " + filename.string());
    }
    for (auto& atom : molecule.atoms) {
        atom.edip_type = atom.type;
        if (atom_type_registry().find(atom.type) == atom_type_registry().end()) {
            register_external_lj(atom.type, 0.0, 0.0);
        }
    }
    for (std::string line; std::getline(input, line);) {
        const auto words = split_ws(strip_comment(line));
        if (words.size() < 18) {
            continue;
        }
        molecule.add_edip_type(words[0], std::stod(words[1]), std::stod(words[2]), std::stod(words[3]),
                               std::stod(words[4]), std::stod(words[5]), std::stod(words[6]),
                               std::stod(words[7]), std::stod(words[8]), std::stod(words[9]),
                               std::stod(words[10]), std::stod(words[11]), std::stod(words[12]),
                               std::stod(words[13]), std::stod(words[14]), std::stod(words[15]),
                               std::stod(words[16]), std::stod(words[17]));
    }
}

std::optional<double> find_external_atom_type_mass(const std::string& atom_type) {
    const auto it = atom_type_registry().find(atom_type);
    if (it == atom_type_registry().end() || it->second.mass <= 0.0) {
        return std::nullopt;
    }
    return it->second.mass;
}

}  // namespace xpongecpp
