#include "core.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace xpongecpp {
namespace {

std::string trim_copy(const std::string& input) {
    const auto first = input.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}

std::vector<std::string> split_ws(const std::string& line) {
    std::istringstream iss(line);
    std::vector<std::string> out;
    std::string word;
    while (iss >> word) {
        out.push_back(word);
    }
    return out;
}

}  // namespace

Assign get_assignment_from_mol2_text(const std::string& text,
                                     std::optional<int> total_charge,
                                     bool total_charge_from_partial_sum) {
    Assign assignment("MOL2");
    std::istringstream input(text);
    std::string line;
    std::string section;
    std::unordered_map<int, std::uint32_t> mol2_to_atom;

    while (std::getline(input, line)) {
        if (line.rfind("@<TRIPOS>", 0) == 0) {
            section = trim_copy(line.substr(9));
            continue;
        }
        if (trim_copy(line).empty()) {
            continue;
        }
        if (section == "MOLECULE" && assignment.name == "MOL2") {
            assignment.name = trim_copy(line);
            continue;
        }
        if (section == "ATOM") {
            const auto words = split_ws(line);
            if (words.size() < 9) {
                throw std::invalid_argument("invalid MOL2 atom line: " + line);
            }
            const int mol2_id = std::stoi(words[0]);
            const auto atom_id = static_cast<std::uint32_t>(assignment.atom_count());
            mol2_to_atom[mol2_id] = atom_id;
            assignment.add_atom(words[5], std::stod(words[2]), std::stod(words[3]), std::stod(words[4]),
                                words[1], std::stod(words[8]));
        } else if (section == "UNITY_ATOM_ATTR") {
            const auto words = split_ws(line);
            if (words.size() < 2) {
                throw std::invalid_argument("invalid MOL2 UNITY_ATOM_ATTR line: " + line);
            }
            const auto atom = mol2_to_atom.at(std::stoi(words[0]));
            const auto count = std::stoi(words[1]);
            for (int i = 0; i < count; ++i) {
                if (!std::getline(input, line)) {
                    throw std::invalid_argument("truncated MOL2 UNITY_ATOM_ATTR section");
                }
                const auto attr_words = split_ws(line);
                if (attr_words.size() < 2) {
                    throw std::invalid_argument("invalid MOL2 UNITY_ATOM_ATTR property line: " + line);
                }
                if (attr_words[0] == "charge") {
                    assignment.formal_charges[atom] = std::stoi(attr_words[1]);
                } else {
                    throw std::invalid_argument("unknown MOL2 UNITY_ATOM_ATTR property: " + attr_words[0]);
                }
            }
        } else if (section == "BOND") {
            const auto words = split_ws(line);
            if (words.size() < 4) {
                throw std::invalid_argument("invalid MOL2 bond line: " + line);
            }
            const auto atom1 = mol2_to_atom.at(std::stoi(words[1]));
            const auto atom2 = mol2_to_atom.at(std::stoi(words[2]));
            if (words[3] == "1" || words[3] == "2" || words[3] == "3" || words[3] == "4" ||
                words[3] == "5" || words[3] == "6" || words[3] == "7" || words[3] == "8" || words[3] == "9") {
                assignment.add_bond(atom1, atom2, std::stoi(words[3]));
            } else if (words[3] == "ar") {
                assignment.add_bond(atom1, atom2, -1);
                assignment.add_bond_marker(atom1, atom2, "mol2_ar");
            } else if (words[3] == "am") {
                assignment.add_bond(atom1, atom2, -1);
                assignment.add_bond_marker(atom1, atom2, "mol2_am");
            } else if (words[3] == "un") {
                assignment.add_bond(atom1, atom2, -1);
            } else {
                throw std::invalid_argument("unsupported MOL2 bond type: " + words[3]);
            }
        }
    }
    if (total_charge_from_partial_sum) {
        double charge_sum = 0.0;
        for (const auto charge : assignment.charges) {
            charge_sum += charge;
        }
        total_charge = static_cast<int>(std::lround(charge_sum));
    }
    const bool success = assignment.determine_bond_order(true, total_charge);
    if (!success) {
        for (auto& bond : assignment.bonds) {
            for (auto& [neighbor, order] : bond) {
                (void)neighbor;
                order = -1;
            }
        }
        assignment.determine_bond_order(true, total_charge);
    }
    return assignment;
}

Assign get_assignment_from_xyz_text(const std::string& text) {
    std::istringstream input(text);
    std::string line;
    std::size_t expected_atoms = 0;
    if (!std::getline(input, line)) {
        throw std::invalid_argument("empty XYZ input");
    }
    try {
        expected_atoms = static_cast<std::size_t>(std::stoul(trim_copy(line)));
    } catch (const std::exception&) {
        throw std::invalid_argument("invalid XYZ atom count: " + line);
    }
    std::getline(input, line);

    Assign assignment("XYZ");
    assignment.elements.reserve(expected_atoms);
    for (std::size_t i = 0; i < expected_atoms; ++i) {
        if (!std::getline(input, line)) {
            throw std::invalid_argument("truncated XYZ input");
        }
        const auto words = split_ws(line);
        if (words.size() < 4) {
            throw std::invalid_argument("invalid XYZ atom line: " + line);
        }
        assignment.add_atom(words[0], std::stod(words[1]), std::stod(words[2]), std::stod(words[3]),
                            words[0] + std::to_string(i + 1));
    }
    assignment.determine_connectivity(1.2);
    return assignment;
}

namespace {

std::string pdb_string(const std::string& line, std::size_t pos, std::size_t len) {
    if (line.size() <= pos) {
        return "";
    }
    return trim_copy(line.substr(pos, std::min(len, line.size() - pos)));
}

double pdb_float(const std::string& line, std::size_t pos, std::size_t len, const char* field_name) {
    if (line.size() <= pos) {
        throw std::invalid_argument(std::string("missing PDB coordinate field: ") + field_name);
    }
    const auto field = trim_copy(line.substr(pos, std::min(len, line.size() - pos)));
    if (field.empty()) {
        throw std::invalid_argument(std::string("empty PDB coordinate field: ") + field_name);
    }
    return std::stod(field);
}

}  // namespace

Assign get_assignment_from_pdb_text(const std::string& text) {
    Assign assignment("PDB");
    std::istringstream input(text);
    std::string line;
    std::unordered_map<int, std::uint32_t> serial_to_atom;
    std::vector<std::pair<int, int>> conect_records;
    while (std::getline(input, line)) {
        if (line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0) {
            const int serial = std::stoi(pdb_string(line, 6, 5));
            const std::string name = pdb_string(line, 12, 4);
            const std::string element = guess_element(name, pdb_string(line, 76, 2));
            const auto atom_id = static_cast<std::uint32_t>(assignment.atom_count());
            serial_to_atom[serial] = atom_id;
            assignment.add_atom(element, pdb_float(line, 30, 8, "x"), pdb_float(line, 38, 8, "y"),
                                pdb_float(line, 46, 8, "z"), name);
        } else if (line.rfind("CONECT", 0) == 0) {
            const auto words = split_ws(line);
            if (words.size() >= 3) {
                const int from = std::stoi(words[1]);
                for (std::size_t i = 2; i < words.size(); ++i) {
                    conect_records.emplace_back(from, std::stoi(words[i]));
                }
            }
        }
    }
    for (const auto& [from, to] : conect_records) {
        const auto from_it = serial_to_atom.find(from);
        const auto to_it = serial_to_atom.find(to);
        if (from_it != serial_to_atom.end() && to_it != serial_to_atom.end() && from_it->second < to_it->second) {
            assignment.add_bond(from_it->second, to_it->second, 1);
        }
    }
    if (assignment.bond_count() == 0 && assignment.atom_count() != 0) {
        assignment.determine_connectivity(1.2);
    }
    return assignment;
}

Assign get_assignment_from_residuetype(const ResidueType& residue_type) {
    Assign assignment(residue_type.name());
    for (const auto& atom : residue_type.atoms()) {
        assignment.add_atom(atom.element.empty() ? guess_element(atom.name, atom.type) : atom.element,
                            atom.x, atom.y, atom.z, atom.name, atom.charge);
    }
    for (const auto& bond : residue_type.bonds()) {
        assignment.add_bond(bond.atom1, bond.atom2, 1);
    }
    return assignment;
}

std::string assignment_to_mol2_text(const Assign& assignment, const std::string& residue_name) {
    std::ostringstream out;
    const std::string resname = residue_name.empty() ? assignment.name : residue_name;
    std::size_t bond_count = 0;
    for (std::uint32_t i = 0; i < assignment.bonds.size(); ++i) {
        for (const auto& [neighbor, order] : assignment.bonds[i]) {
            (void)order;
            if (neighbor > i) {
                ++bond_count;
            }
        }
    }
    out << "@<TRIPOS>MOLECULE\n" << assignment.name << "\n";
    out << std::setw(6) << assignment.atom_count() << std::setw(6) << bond_count << std::setw(6) << 1
        << "     0     1\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n";
    out << std::fixed << std::setprecision(4);
    for (std::uint32_t i = 0; i < assignment.atom_count(); ++i) {
        const auto& coord = assignment.coordinates[i];
        const std::string atom_name = assignment.names[i].empty() ? assignment.elements[i] + std::to_string(i + 1)
                                                                  : assignment.names[i];
        const std::string atom_type = i < assignment.atom_types.size() && !assignment.atom_types[i].empty()
                                          ? assignment.atom_types[i]
                                          : assignment.elements[i] + assignment.element_details[i];
        out << std::setw(6) << i + 1 << " " << std::setw(4) << atom_name << " "
            << std::setw(10) << coord[0] << " " << std::setw(10) << coord[1] << " " << std::setw(10) << coord[2]
            << " " << std::setw(6) << atom_type << " " << std::setw(5) << 1 << " " << std::setw(8) << resname
            << " " << std::setw(10) << std::setprecision(6) << assignment.charges[i] << std::setprecision(4)
            << "\n";
    }
    out << "@<TRIPOS>BOND\n";
    std::size_t bond_index = 1;
    for (std::uint32_t i = 0; i < assignment.bonds.size(); ++i) {
        for (const auto& [neighbor, order] : assignment.bonds[i]) {
            if (neighbor <= i) {
                continue;
            }
            out << std::setw(6) << bond_index++ << std::setw(6) << i + 1 << std::setw(6) << neighbor + 1
                << " " << (order > 0 ? order : 1) << "\n";
        }
    }
    out << "@<TRIPOS>SUBSTRUCTURE\n" << std::setw(5) << 1 << " " << std::setw(8) << resname << std::setw(6)
        << 1 << " ****               0 ****  **** \n";
    return out.str();
}

std::string assignment_to_pdb_text(const Assign& assignment, const std::string& residue_name) {
    std::ostringstream out;
    const std::string resname = residue_name.empty() ? assignment.name : residue_name;
    out << std::fixed << std::setprecision(3);
    for (std::uint32_t i = 0; i < assignment.atom_count(); ++i) {
        const auto& coord = assignment.coordinates[i];
        const std::string atom_name = assignment.names[i].empty() ? assignment.elements[i] + std::to_string(i + 1)
                                                                  : assignment.names[i];
        out << "ATOM  " << std::setw(5) << i + 1 << " " << std::left << std::setw(4) << atom_name
            << std::right << " " << std::setw(3) << resname << " A" << std::setw(4) << 1 << "    "
            << std::setw(8) << coord[0] << std::setw(8) << coord[1] << std::setw(8) << coord[2]
            << "  1.00  0.00          " << std::setw(2) << assignment.elements[i] << "\n";
    }
    out << "END\n";
    return out.str();
}

}  // namespace xpongecpp
