#include "core.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <stdexcept>

namespace xpongecpp {

Assign::Assign(std::string assign_name) : name(std::move(assign_name)) {}

void Assign::add_atom(const std::string& element, double x, double y, double z,
                      const std::string& atom_name, double charge) {
    std::string base_element = element;
    std::string detail;
    const auto dot = element.find('.');
    if (dot != std::string::npos) {
        base_element = element.substr(0, dot);
        detail = element.substr(dot);
    }
    elements.push_back(base_element);
    element_details.push_back(detail);
    names.push_back(atom_name.empty() ? element : atom_name);
    coordinates.push_back({x, y, z});
    charges.push_back(charge);
    formal_charges.push_back(0);
    bonds.emplace_back();
    bond_markers.emplace_back();
    atom_markers.emplace_back();
    atom_types.emplace_back();
    built = false;
}

void Assign::add_bond(std::uint32_t atom1, std::uint32_t atom2, int order) {
    if (atom1 >= elements.size() || atom2 >= elements.size()) {
        throw std::out_of_range("Assign bond atom index out of range");
    }
    if (atom1 == atom2) {
        throw std::invalid_argument("Assign self bond");
    }
    bonds[atom1][atom2] = order;
    bonds[atom2][atom1] = order;
    bond_markers[atom1][atom2] = {};
    bond_markers[atom2][atom1] = {};
    built = false;
}

void Assign::set_charge(std::uint32_t atom, double charge) {
    if (atom >= charges.size()) {
        throw std::out_of_range("Assign charge atom index out of range");
    }
    charges[atom] = charge;
}

void Assign::set_charges(const std::vector<double>& new_charges) {
    if (new_charges.size() != charges.size()) {
        throw std::invalid_argument("Assign charge count does not match atom count");
    }
    charges = new_charges;
}

void Assign::set_formal_charge(std::uint32_t atom, int charge) {
    if (atom >= formal_charges.size()) {
        throw std::out_of_range("Assign formal charge atom index out of range");
    }
    formal_charges[atom] = charge;
    built = false;
}

void Assign::determine_connectivity(double simple_cutoff) {
    if (simple_cutoff <= 0.0) {
        throw std::invalid_argument("simple_cutoff should be positive");
    }
    const double cutoff2 = simple_cutoff * simple_cutoff;
    for (std::uint32_t i = 0; i < coordinates.size(); ++i) {
        for (std::uint32_t j = i + 1; j < coordinates.size(); ++j) {
            const double dx = coordinates[i][0] - coordinates[j][0];
            const double dy = coordinates[i][1] - coordinates[j][1];
            const double dz = coordinates[i][2] - coordinates[j][2];
            const double d2 = dx * dx + dy * dy + dz * dz;
            if (d2 < cutoff2) {
                add_bond(i, j, 1);
            }
        }
    }
}

bool Assign::check_connectivity() const {
    if (elements.empty()) {
        return true;
    }
    std::vector<std::uint8_t> visited(elements.size(), 0);
    std::vector<std::uint32_t> stack{0};
    visited[0] = 1;
    while (!stack.empty()) {
        const auto atom = stack.back();
        stack.pop_back();
        for (const auto& [neighbor, order] : bonds[atom]) {
            (void)order;
            if (visited[neighbor] == 0) {
                visited[neighbor] = 1;
                stack.push_back(neighbor);
            }
        }
    }
    return std::all_of(visited.begin(), visited.end(), [](std::uint8_t value) { return value != 0; });
}

bool Assign::atom_judge(std::uint32_t atom, const std::string& mask) const {
    if (atom >= elements.size()) {
        return false;
    }
    std::string element;
    std::string digits;
    for (const char c : mask) {
        if (std::isdigit(static_cast<unsigned char>(c))) {
            digits.push_back(c);
        } else {
            element.push_back(c);
        }
    }
    if (digits.empty()) {
        return elements[atom] == element;
    }
    return elements[atom] == element && bonds[atom].size() == static_cast<std::size_t>(std::stoi(digits));
}

bool Assign::atom_judge(std::uint32_t atom, const std::vector<std::string>& masks) const {
    return std::any_of(masks.begin(), masks.end(), [&](const std::string& mask) { return atom_judge(atom, mask); });
}

void Assign::add_atom_marker(std::uint32_t atom, const std::string& marker) {
    if (atom >= atom_markers.size()) {
        throw std::out_of_range("Assign atom marker atom index out of range");
    }
    ++atom_markers[atom][marker];
}

void Assign::add_bond_marker(std::uint32_t atom1, std::uint32_t atom2, const std::string& marker, bool only1) {
    if (atom1 >= bond_markers.size() || atom2 >= bond_markers.size()) {
        throw std::out_of_range("Assign bond marker atom index out of range");
    }
    bond_markers[atom1][atom2].insert(marker);
    add_atom_marker(atom1, marker);
    if (!only1) {
        bond_markers[atom2][atom1].insert(marker);
        add_atom_marker(atom2, marker);
    }
}

bool Assign::has_atom_marker(std::uint32_t atom, const std::string& marker) const {
    if (atom >= atom_markers.size()) {
        return false;
    }
    return atom_markers[atom].count(marker) != 0;
}

int Assign::atom_marker_count(std::uint32_t atom, const std::string& marker) const {
    if (atom >= atom_markers.size()) {
        return 0;
    }
    const auto found = atom_markers[atom].find(marker);
    return found == atom_markers[atom].end() ? 0 : found->second;
}

bool Assign::has_bond_marker(std::uint32_t atom1, std::uint32_t atom2, const std::string& marker) const {
    if (atom1 >= bond_markers.size()) {
        return false;
    }
    const auto found = bond_markers[atom1].find(atom2);
    if (found == bond_markers[atom1].end()) {
        return false;
    }
    return found->second.count(marker) != 0;
}

ResidueType Assign::to_residuetype(const std::string& residue_name) const {
    ResidueType residue_type(residue_name);
    for (std::size_t i = 0; i < elements.size(); ++i) {
        const auto& atom_type = atom_types[i].empty() ? elements[i] : atom_types[i];
        residue_type.add_atom(names[i], atom_type, coordinates[i][0], coordinates[i][1], coordinates[i][2],
                              charges[i], default_mass_for_element(elements[i]));
    }
    for (std::size_t i = 0; i < bonds.size(); ++i) {
        for (const auto& [j, order] : bonds[i]) {
            (void)order;
            if (i < j) {
                residue_type.add_connectivity(names[i], names[j]);
            }
        }
    }
    return residue_type;
}

Molecule Assign::to_molecule(const std::string& residue_name) const {
    Molecule molecule(residue_name);
    molecule.append_residue_from_type(to_residuetype(residue_name), 0.0, 0.0, 0.0);
    return molecule;
}

std::size_t Assign::atom_count() const noexcept {
    return elements.size();
}

std::size_t Assign::bond_count() const noexcept {
    std::size_t count = 0;
    for (const auto& atom_bonds : bonds) {
        count += atom_bonds.size();
    }
    return count / 2;
}

}  // namespace xpongecpp
