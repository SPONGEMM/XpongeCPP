#include "nonamber_internal.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace xpongecpp {

std::unordered_map<std::string, AtomTypeInfo>& atom_type_registry() {
    static std::unordered_map<std::string, AtomTypeInfo> registry;
    return registry;
}

std::vector<CharmmUreyParameter>& charmm_urey_registry() {
    static std::vector<CharmmUreyParameter> registry;
    return registry;
}

std::vector<NBfixParameter>& nbfix_registry() {
    static std::vector<NBfixParameter> registry;
    return registry;
}

std::string trim_copy(const std::string& input) {
    const auto first = input.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}

std::string strip_comment(const std::string& line) {
    const auto semicolon = line.find(';');
    const auto hash = line.find('#');
    auto cut = std::min(semicolon == std::string::npos ? line.size() : semicolon,
                        hash == std::string::npos ? line.size() : hash);
    if (line.rfind("#include", 0) == 0) {
        cut = line.size();
    }
    return trim_copy(line.substr(0, cut));
}

std::vector<std::string> split_ws(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> words;
    std::string word;
    while (stream >> word) {
        words.push_back(word);
    }
    return words;
}

std::string upper_copy(std::string text) {
    std::transform(text.begin(), text.end(), text.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
    });
    return text;
}

AtomId one_based_atom(const std::string& token) {
    const auto value = std::stoul(token);
    if (value == 0) {
        throw std::out_of_range("atom indices are one-based in external topology files");
    }
    return static_cast<AtomId>(value - 1);
}

void register_external_lj(const std::string& type, double epsilon, double rmin) {
    atom_type_registry()[type].epsilon = epsilon;
    atom_type_registry()[type].rmin = rmin;
    register_amber_lj_parameter(type, type, epsilon, rmin);
}

void register_external_mass(const std::string& type, double mass) {
    atom_type_registry()[type].mass = mass;
}

void append_atom(Molecule& molecule, std::unordered_map<int, ResidueId>& residue_id_by_number,
                 int residue_number, const std::string& residue_name, const std::string& atom_name,
                 const std::string& atom_type, double charge, double mass) {
    auto residue_it = residue_id_by_number.find(residue_number);
    if (residue_it == residue_id_by_number.end()) {
        const auto residue_id = static_cast<ResidueId>(molecule.residues.size());
        residue_it = residue_id_by_number.emplace(residue_number, residue_id).first;
        Residue residue;
        residue.name = residue_name;
        residue.type_name = residue_name;
        residue.original_name = residue_name;
        residue.atom_begin = static_cast<AtomId>(molecule.atoms.size());
        molecule.residues.push_back(std::move(residue));
    }
    auto& residue = molecule.residues[residue_it->second];
    Atom atom;
    atom.name = atom_name;
    atom.type = atom_type;
    atom.element = mass > 0.0 ? guess_element_from_mass(mass) : guess_element(atom_name, "");
    atom.residue = residue_it->second;
    atom.charge = charge;
    atom.mass = mass;
    molecule.atoms.push_back(std::move(atom));
    ++residue.atom_count;
}

namespace {

void parse_gromacs_file(const std::filesystem::path& filename, Molecule& molecule,
                        std::unordered_set<std::filesystem::path>& included) {
    const auto canonical = std::filesystem::absolute(filename);
    if (!included.insert(canonical).second) {
        return;
    }
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open GROMACS topology file: " + filename.string());
    }

    std::string section;
    std::unordered_map<int, ResidueId> residue_id_by_number;
    for (std::string line; std::getline(input, line);) {
        const auto trimmed = trim_copy(line);
        if (trimmed.rfind("#include", 0) == 0) {
            const auto first_quote = trimmed.find('"');
            const auto last_quote = trimmed.find_last_of('"');
            if (first_quote == std::string::npos || last_quote <= first_quote) {
                throw std::runtime_error("malformed #include line: " + trimmed);
            }
            parse_gromacs_file(filename.parent_path() / trimmed.substr(first_quote + 1, last_quote - first_quote - 1),
                               molecule, included);
            continue;
        }

        const auto clean = strip_comment(line);
        if (clean.empty()) {
            continue;
        }
        if (clean.front() == '[' && clean.back() == ']') {
            section = upper_copy(trim_copy(clean.substr(1, clean.size() - 2)));
            static const std::unordered_set<std::string> supported{
                "DEFAULTS", "ATOMTYPES", "BONDTYPES", "CONSTRAINTTYPES", "ANGLETYPES", "DIHEDRALTYPES",
                "PAIRTYPES", "NONBOND_PARAMS", "MOLECULETYPE", "ATOMS", "BONDS", "PAIRS", "ANGLES",
                "DIHEDRALS", "CMAPTYPES", "CMAP", "VIRTUAL_SITES3", "EXCLUSIONS", "SYSTEM", "MOLECULES"};
            if (supported.find(section) == supported.end()) {
                throw std::runtime_error("unsupported GROMACS topology section: " + section);
            }
            continue;
        }

        const auto words = split_ws(clean);
        if (words.empty()) {
            continue;
        }
        if (section == "ATOMTYPES") {
            if (words.size() < 7) {
                throw std::runtime_error("malformed atomtypes row");
            }
            const auto type = words[0];
            const auto mass = words[4] == "A" ? std::stod(words[2]) : std::stod(words[3]);
            const auto sigma = std::stod(words[words.size() - 2]);
            const auto epsilon = std::stod(words[words.size() - 1]);
            register_external_mass(type, mass);
            register_external_lj(type, epsilon, sigma / 2.0);
        } else if (section == "PAIRTYPES" || section == "NONBOND_PARAMS") {
            if (words.size() >= 5) {
                nbfix_registry().push_back({words[0], words[1], std::stod(words[3]), std::stod(words[4]), 0.0});
            }
        } else if (section == "ANGLETYPES") {
            if (words.size() >= 8 && words[3] == "5") {
                charmm_urey_registry().push_back({{words[0], words[1], words[2]},
                                                 std::stod(words[5]), std::stod(words[4]),
                                                 std::stod(words[7]), std::stod(words[6])});
            }
        } else if (section == "MOLECULETYPE") {
            molecule.name = words[0];
        } else if (section == "ATOMS") {
            if (words.size() < 8) {
                throw std::runtime_error("malformed atoms row");
            }
            append_atom(molecule, residue_id_by_number, std::stoi(words[2]), words[3], words[4], words[1],
                        std::stod(words[6]), std::stod(words[7]));
        } else if (section == "BONDS") {
            if (words.size() < 2) {
                throw std::runtime_error("malformed bonds row");
            }
            molecule.explicit_bonds.push_back({one_based_atom(words[0]), one_based_atom(words[1])});
        } else if (section == "PAIRS") {
            if (words.size() >= 6) {
                molecule.add_nb14_extra(one_based_atom(words[0]), one_based_atom(words[1]),
                                        std::stod(words[3]), std::stod(words[4]), std::stod(words[5]));
            }
        } else if (section == "ANGLES") {
            if (words.size() >= 8 && words[3] == "5") {
                molecule.add_urey_bradley(one_based_atom(words[0]), one_based_atom(words[1]), one_based_atom(words[2]),
                                          std::stod(words[5]), std::stod(words[4]), std::stod(words[7]),
                                          std::stod(words[6]));
            }
        } else if (section == "DIHEDRALS") {
            if (words.size() >= 11 && words[4] == "3") {
                molecule.add_ryckaert_bellemans(one_based_atom(words[0]), one_based_atom(words[1]),
                                                one_based_atom(words[2]), one_based_atom(words[3]),
                                                std::stod(words[5]), std::stod(words[6]), std::stod(words[7]),
                                                std::stod(words[8]), std::stod(words[9]), std::stod(words[10]));
            }
        } else if (section == "VIRTUAL_SITES3") {
            if (words.size() >= 7) {
                molecule.add_virtual_atom2(one_based_atom(words[0]), one_based_atom(words[1]),
                                           one_based_atom(words[2]), one_based_atom(words[3]),
                                           std::stod(words[5]), std::stod(words[6]));
            }
        } else if (section == "CMAP") {
            if (words.size() >= 8) {
                const auto resolution = static_cast<std::uint32_t>(std::stoul(words[6]));
                std::vector<double> parameters;
                parameters.reserve(resolution * resolution);
                for (std::size_t i = 7; i < words.size(); ++i) {
                    parameters.push_back(std::stod(words[i]));
                }
                const auto type = molecule.add_cmap_type(resolution, parameters);
                molecule.add_cmap(one_based_atom(words[0]), one_based_atom(words[1]), one_based_atom(words[2]),
                                  one_based_atom(words[3]), one_based_atom(words[4]), type);
            }
        }
    }
}

}  // namespace

Molecule load_gromacs_topology_file(const std::filesystem::path& filename) {
    Molecule molecule("GMX");
    std::unordered_set<std::filesystem::path> included;
    parse_gromacs_file(filename, molecule, included);
    return molecule;
}

}  // namespace xpongecpp
