#include "core.hpp"

#include <cmath>
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

}  // namespace xpongecpp
