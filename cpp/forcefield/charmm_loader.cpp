#include "nonamber_internal.hpp"

#include <fstream>
#include <optional>
#include <stdexcept>

namespace xpongecpp {
namespace {

void add_nbfixes_for_atoms(Molecule& molecule) {
    for (AtomId atom1 = 0; atom1 < molecule.atoms.size(); ++atom1) {
        for (AtomId atom2 = atom1 + 1; atom2 < molecule.atoms.size(); ++atom2) {
            for (const auto& nbfix : nbfix_registry()) {
                const auto& type1 = molecule.atoms[atom1].type;
                const auto& type2 = molecule.atoms[atom2].type;
                if ((type1 == nbfix.type1 && type2 == nbfix.type2) ||
                    (type1 == nbfix.type2 && type2 == nbfix.type1)) {
                    molecule.add_nb14_extra(atom1, atom2, nbfix.a, nbfix.b, nbfix.kee);
                }
            }
        }
    }
}

std::optional<CharmmUreyParameter> find_charmm_urey(const std::string& type1, const std::string& type2,
                                                    const std::string& type3) {
    for (const auto& parameter : charmm_urey_registry()) {
        if ((parameter.types[0] == type1 && parameter.types[1] == type2 && parameter.types[2] == type3) ||
            (parameter.types[0] == type3 && parameter.types[1] == type2 && parameter.types[2] == type1)) {
            return parameter;
        }
    }
    return std::nullopt;
}

std::vector<std::vector<AtomId>> adjacency_from_explicit_bonds(const Molecule& molecule) {
    std::vector<std::vector<AtomId>> adjacency(molecule.atoms.size());
    for (const auto& bond : molecule.explicit_bonds) {
        adjacency[bond.atom1].push_back(bond.atom2);
        adjacency[bond.atom2].push_back(bond.atom1);
    }
    return adjacency;
}

void add_charmm_ureys_for_angles(Molecule& molecule) {
    const auto adjacency = adjacency_from_explicit_bonds(molecule);
    for (AtomId center = 0; center < adjacency.size(); ++center) {
        const auto& neighbors = adjacency[center];
        for (std::size_t i = 0; i < neighbors.size(); ++i) {
            for (std::size_t k = i + 1; k < neighbors.size(); ++k) {
                const auto match = find_charmm_urey(molecule.atoms[neighbors[i]].type,
                                                    molecule.atoms[center].type,
                                                    molecule.atoms[neighbors[k]].type);
                if (match) {
                    molecule.add_urey_bradley(neighbors[i], center, neighbors[k],
                                              match->k, match->theta, match->k_ub, match->r13);
                }
            }
        }
    }
}

}  // namespace

void load_charmm_parameter_file(const std::filesystem::path& filename) {
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open CHARMM parameter file: " + filename.string());
    }
    std::string section;
    for (std::string line; std::getline(input, line);) {
        const auto clean = strip_comment(line);
        if (clean.empty()) {
            continue;
        }
        const auto words = split_ws(clean);
        if (words.empty()) {
            continue;
        }
        const auto first = upper_copy(words[0]);
        if (first == "MASS" && words.size() >= 4) {
            register_external_mass(words[2], std::stod(words[3]));
            continue;
        }
        if (first == "NONBONDED" || first == "ANGLE" || first == "ANGL" || first == "NBFIX") {
            section = first == "ANGL" ? "ANGLE" : first;
            continue;
        }
        if (section == "NONBONDED" && words.size() >= 4) {
            register_external_lj(words[0], std::abs(std::stod(words[2])), std::stod(words[3]));
        } else if (section == "ANGLE" && words.size() >= 7) {
            charmm_urey_registry().push_back({{words[0], words[1], words[2]},
                                             std::stod(words[3]), std::stod(words[4]),
                                             std::stod(words[5]), std::stod(words[6])});
        } else if (section == "NBFIX" && words.size() >= 5) {
            const double kee = words.size() >= 6 ? std::stod(words[5]) : 0.0;
            nbfix_registry().push_back({words[0], words[1], std::stod(words[2]), std::stod(words[3]), kee});
        }
    }
}

Molecule load_charmm_topology_file(const std::filesystem::path& filename) {
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open CHARMM topology file: " + filename.string());
    }
    Molecule molecule("CHARMM");
    std::unordered_map<std::string, AtomId> atom_by_name;
    std::unordered_map<int, ResidueId> residue_id_by_number;
    int residue_number = 0;
    std::string residue_name = "RES";
    for (std::string line; std::getline(input, line);) {
        const auto clean = strip_comment(line);
        if (clean.empty()) {
            continue;
        }
        const auto words = split_ws(clean);
        const auto key = upper_copy(words[0]);
        if (key == "RESI" || key == "PRES") {
            residue_name = words.at(1);
            residue_number = static_cast<int>(molecule.residues.size()) + 1;
            atom_by_name.clear();
        } else if (key == "ATOM" && words.size() >= 4) {
            const auto type = words[2];
            const auto mass_it = atom_type_registry().find(type);
            const double mass = mass_it == atom_type_registry().end() ? 0.0 : mass_it->second.mass;
            append_atom(molecule, residue_id_by_number, residue_number, residue_name, words[1], type,
                        std::stod(words[3]), mass);
            atom_by_name[words[1]] = static_cast<AtomId>(molecule.atoms.size() - 1);
        } else if (key == "BOND") {
            for (std::size_t i = 1; i + 1 < words.size(); i += 2) {
                molecule.explicit_bonds.push_back({atom_by_name.at(words[i]), atom_by_name.at(words[i + 1])});
            }
        }
    }
    add_charmm_ureys_for_angles(molecule);
    add_nbfixes_for_atoms(molecule);
    return molecule;
}

}  // namespace xpongecpp
