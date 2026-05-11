#include "core.hpp"

#include <algorithm>
#include <sstream>
#include <stdexcept>

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

Molecule load_mol2_text(const std::string& text) {
    Molecule molecule("MOL2");
    std::istringstream input(text);
    std::string line;
    std::string section;
    std::unordered_map<int, AtomId> mol2_to_atom;
    std::unordered_map<int, ResidueId> residue_index_to_id;

    while (std::getline(input, line)) {
        if (line.rfind("@<TRIPOS>", 0) == 0) {
            section = trim_copy(line.substr(9));
            continue;
        }
        if (trim_copy(line).empty()) {
            continue;
        }
        if (section == "MOLECULE" && molecule.name == "MOL2") {
            molecule.name = trim_copy(line);
            continue;
        }
        if (section == "ATOM") {
            const auto words = split_ws(line);
            if (words.size() < 9) {
                throw std::invalid_argument("invalid MOL2 atom line: " + line);
            }
            const int mol2_id = std::stoi(words[0]);
            const int residue_index = std::stoi(words[6]);
            const std::string residue_name = words[7];
            ResidueId residue_id = 0;
            const auto residue_it = residue_index_to_id.find(residue_index);
            if (residue_it == residue_index_to_id.end()) {
                residue_id = static_cast<ResidueId>(molecule.residues.size());
                residue_index_to_id[residue_index] = residue_id;
                Residue residue;
                residue.name = residue_name;
                residue.type_name = residue_name;
                residue.atom_begin = static_cast<AtomId>(molecule.atoms.size());
                residue.atom_count = 0;
                molecule.residues.push_back(residue);
            } else {
                residue_id = residue_it->second;
            }

            Atom atom;
            atom.name = words[1];
            atom.x = std::stod(words[2]);
            atom.y = std::stod(words[3]);
            atom.z = std::stod(words[4]);
            atom.type = words[5];
            atom.element = guess_element(atom.name, atom.type);
            atom.residue = residue_id;
            atom.charge = std::stod(words[8]);
            atom.mass = default_mass_for_element(atom.element);
            mol2_to_atom[mol2_id] = static_cast<AtomId>(molecule.atoms.size());
            molecule.atoms.push_back(std::move(atom));
            molecule.residues[residue_id].atom_count += 1;
        } else if (section == "BOND") {
            const auto words = split_ws(line);
            if (words.size() < 4) {
                throw std::invalid_argument("invalid MOL2 bond line: " + line);
            }
            const auto atom1 = mol2_to_atom.at(std::stoi(words[1]));
            const auto atom2 = mol2_to_atom.at(std::stoi(words[2]));
            molecule.explicit_bonds.push_back({atom1, atom2});
            if (molecule.atoms[atom1].residue != molecule.atoms[atom2].residue) {
                molecule.residue_links.push_back({atom1, atom2});
            }
        }
    }
    return molecule;
}

}  // namespace xpongecpp
