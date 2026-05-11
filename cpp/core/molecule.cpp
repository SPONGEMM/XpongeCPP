#include "core.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
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

void ensure_atom_id(const Molecule& molecule, AtomId id) {
    if (id >= molecule.atoms.size()) {
        throw std::out_of_range("AtomId out of range");
    }
}

void ensure_residue_id(const Molecule& molecule, ResidueId id) {
    if (id >= molecule.residues.size()) {
        throw std::out_of_range("ResidueId out of range");
    }
}

std::array<double, 3> molecule_min(const Molecule& molecule) {
    std::array<double, 3> minv{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
    };
    for (const auto& atom : molecule.atoms) {
        minv[0] = std::min(minv[0], atom.x);
        minv[1] = std::min(minv[1], atom.y);
        minv[2] = std::min(minv[2], atom.z);
    }
    return minv;
}

std::array<double, 3> molecule_max(const Molecule& molecule) {
    std::array<double, 3> maxv{
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
    };
    for (const auto& atom : molecule.atoms) {
        maxv[0] = std::max(maxv[0], atom.x);
        maxv[1] = std::max(maxv[1], atom.y);
        maxv[2] = std::max(maxv[2], atom.z);
    }
    return maxv;
}

}  // namespace

double default_mass_for_element(const std::string& element) {
    if (element == "H") return 1.008;
    if (element == "C") return 12.010;
    if (element == "N") return 14.010;
    if (element == "O") return 16.000;
    if (element == "P") return 30.970;
    if (element == "S") return 32.060;
    if (element == "Cl" || element == "CL") return 35.450;
    if (element == "Na" || element == "NA") return 22.990;
    if (element == "K") return 39.100;
    return 0.0;
}

std::string guess_element(const std::string& atom_name, const std::string& explicit_element) {
    const auto explicit_trimmed = trim_copy(explicit_element);
    if (!explicit_trimmed.empty()) {
        if (!std::isalpha(static_cast<unsigned char>(explicit_trimmed[0]))) {
            return guess_element(atom_name, "");
        }
        if (explicit_trimmed == "OW" || explicit_trimmed == "OH" || explicit_trimmed == "O2") return "O";
        if (explicit_trimmed == "HW" || explicit_trimmed == "HO" || explicit_trimmed[0] == 'H') return "H";
        if (explicit_trimmed == "Na+" || explicit_trimmed == "NA") return "Na";
        if (explicit_trimmed == "Cl-" || explicit_trimmed == "CL") return "Cl";
        if (std::islower(static_cast<unsigned char>(explicit_trimmed[0]))) {
            if (explicit_trimmed.rfind("cl", 0) == 0) return "Cl";
            if (explicit_trimmed.rfind("br", 0) == 0) return "Br";
            return std::string(1, static_cast<char>(std::toupper(static_cast<unsigned char>(explicit_trimmed[0]))));
        }
        if (explicit_trimmed.size() >= 2 && std::islower(static_cast<unsigned char>(explicit_trimmed[1]))) {
            std::string element = explicit_trimmed.substr(0, 2);
            element[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(element[0])));
            return element;
        }
        return std::string(1, static_cast<char>(std::toupper(static_cast<unsigned char>(explicit_trimmed[0]))));
    }
    std::string letters;
    for (const char c : atom_name) {
        if (std::isalpha(static_cast<unsigned char>(c))) {
            letters.push_back(c);
        }
    }
    if (letters.empty()) {
        return "X";
    }
    std::string element;
    element.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(letters[0]))));
    if (letters.size() > 1 && std::islower(static_cast<unsigned char>(letters[1]))) {
        element.push_back(letters[1]);
    }
    return element;
}

ResidueType::ResidueType(std::string name) : name_(std::move(name)) {}

const std::string& ResidueType::name() const noexcept { return name_; }
std::uint64_t ResidueType::version() const noexcept { return version_; }
std::size_t ResidueType::atom_count() const noexcept { return atoms_.size(); }
std::size_t ResidueType::bond_count() const noexcept { return bonds_.size(); }
const std::vector<ResidueTypeAtom>& ResidueType::atoms() const noexcept { return atoms_; }
const std::vector<ResidueTypeBond>& ResidueType::bonds() const noexcept { return bonds_; }

void ResidueType::add_atom(const std::string& name, const std::string& type, double x, double y, double z,
                           double charge, double mass) {
    if (atom_name_to_index_.count(name) != 0) {
        throw std::invalid_argument("duplicate atom name in ResidueType: " + name);
    }
    ResidueTypeAtom atom;
    atom.name = name;
    atom.type = type;
    atom.element = guess_element(name, type);
    atom.x = x;
    atom.y = y;
    atom.z = z;
    atom.charge = charge;
    atom.mass = mass == 0.0 ? default_mass_for_element(atom.element) : mass;
    atom_name_to_index_[name] = static_cast<std::uint32_t>(atoms_.size());
    atoms_.push_back(std::move(atom));
    ++version_;
}

void ResidueType::add_connectivity(const std::string& atom1, const std::string& atom2) {
    const auto index1 = atom_index(atom1);
    const auto index2 = atom_index(atom2);
    if (index1 == index2) {
        throw std::invalid_argument("self bond in ResidueType");
    }
    const auto lo = std::min(index1, index2);
    const auto hi = std::max(index1, index2);
    for (const auto& bond : bonds_) {
        if (bond.atom1 == lo && bond.atom2 == hi) {
            return;
        }
    }
    bonds_.push_back({lo, hi});
    ++version_;
}

std::uint32_t ResidueType::atom_index(const std::string& name) const {
    const auto it = atom_name_to_index_.find(name);
    if (it == atom_name_to_index_.end()) {
        throw std::out_of_range("atom name not found in ResidueType: " + name);
    }
    return it->second;
}

Molecule::Molecule(std::string molecule_name) : name(std::move(molecule_name)) {}
std::size_t Molecule::atom_count() const noexcept { return atoms.size(); }
std::size_t Molecule::residue_count() const noexcept { return residues.size(); }

const Atom& Molecule::atom(AtomId id) const {
    ensure_atom_id(*this, id);
    return atoms[id];
}

Atom& Molecule::atom(AtomId id) {
    ensure_atom_id(*this, id);
    return atoms[id];
}

const Residue& Molecule::residue(ResidueId id) const {
    ensure_residue_id(*this, id);
    return residues[id];
}

Residue& Molecule::residue(ResidueId id) {
    ensure_residue_id(*this, id);
    return residues[id];
}

void Molecule::append_residue_from_type(const ResidueType& type, double dx, double dy, double dz) {
    const ResidueId residue_id = static_cast<ResidueId>(residues.size());
    Residue residue;
    residue.name = type.name();
    residue.type_name = type.name();
    residue.atom_begin = static_cast<AtomId>(atoms.size());
    residue.atom_count = static_cast<std::uint32_t>(type.atom_count());
    residues.push_back(residue);

    for (const auto& template_atom : type.atoms()) {
        Atom atom;
        atom.name = template_atom.name;
        atom.type = template_atom.type;
        atom.element = template_atom.element;
        atom.residue = residue_id;
        atom.x = template_atom.x + dx;
        atom.y = template_atom.y + dy;
        atom.z = template_atom.z + dz;
        atom.charge = template_atom.charge;
        atom.mass = template_atom.mass;
        atoms.push_back(std::move(atom));
    }
    for (const auto& bond : type.bonds()) {
        explicit_bonds.push_back({residue.atom_begin + bond.atom1, residue.atom_begin + bond.atom2});
    }
}

void Molecule::set_box_padding(double padding, bool center) {
    if (padding < 0.0) {
        throw std::invalid_argument("padding should be non-negative");
    }
    if (atoms.empty()) {
        throw std::invalid_argument("at least one atom is required to set box padding");
    }
    const auto minv = molecule_min(*this);
    const auto maxv = molecule_max(*this);
    box_length = {
        maxv[0] - minv[0] + 2.0 * padding,
        maxv[1] - minv[1] + 2.0 * padding,
        maxv[2] - minv[2] + 2.0 * padding,
    };
    has_box = true;
    if (center) {
        const std::array<double, 3> shift{padding - minv[0], padding - minv[1], padding - minv[2]};
        for (auto& atom : atoms) {
            atom.x += shift[0];
            atom.y += shift[1];
            atom.z += shift[2];
        }
    }
}

bool Molecule::validate() const {
    for (std::size_t residue_index = 0; residue_index < residues.size(); ++residue_index) {
        const auto& res = residues[residue_index];
        if (static_cast<std::size_t>(res.atom_begin) + res.atom_count > atoms.size()) {
            return false;
        }
        for (std::uint32_t local = 0; local < res.atom_count; ++local) {
            const auto atom_id = res.atom_begin + local;
            if (atoms[atom_id].residue != residue_index) {
                return false;
            }
        }
    }
    for (const auto& link : residue_links) {
        if (link.atom1 >= atoms.size() || link.atom2 >= atoms.size()) {
            return false;
        }
    }
    for (const auto& bond : explicit_bonds) {
        if (bond.atom1 >= atoms.size() || bond.atom2 >= atoms.size() || bond.atom1 == bond.atom2) {
            return false;
        }
    }
    return true;
}

std::unordered_map<std::string, std::size_t> Molecule::residue_counts() const {
    std::unordered_map<std::string, std::size_t> counts;
    for (const auto& residue : residues) {
        counts[residue.name] += 1;
    }
    return counts;
}

Assign::Assign(std::string assign_name) : name(std::move(assign_name)) {}

void Assign::add_atom(const std::string& element, double x, double y, double z,
                      const std::string& atom_name, double charge) {
    elements.push_back(element);
    names.push_back(atom_name.empty() ? element : atom_name);
    coordinates.push_back({x, y, z});
    charges.push_back(charge);
    bonds.emplace_back();
    atom_types.emplace_back();
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

void Assign::determine_atom_type(const std::string&) {
    for (std::size_t i = 0; i < elements.size(); ++i) {
        atom_types[i] = elements[i];
    }
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

}  // namespace xpongecpp
