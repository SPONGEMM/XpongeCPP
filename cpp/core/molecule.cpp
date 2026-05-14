#include "core.hpp"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <limits>
#include <stdexcept>

namespace xpongecpp {
namespace {

constexpr AtomId invalid_atom_id = std::numeric_limits<AtomId>::max();

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

AtomId remap_atom_id(const std::vector<AtomId>& old_to_new_atom, AtomId old_atom_id) {
    if (old_atom_id >= old_to_new_atom.size()) {
        return invalid_atom_id;
    }
    return old_to_new_atom[old_atom_id];
}

void append_internal_structures(Molecule& target, const Molecule& source, AtomId atom_offset) {
    target.explicit_bonds.reserve(target.explicit_bonds.size() + source.explicit_bonds.size());
    for (const auto& bond : source.explicit_bonds) {
        target.explicit_bonds.push_back({bond.atom1 + atom_offset, bond.atom2 + atom_offset});
    }
    target.residue_links.reserve(target.residue_links.size() + source.residue_links.size());
    for (const auto& link : source.residue_links) {
        target.residue_links.push_back({link.atom1 + atom_offset, link.atom2 + atom_offset});
    }
    target.virtual_atoms.reserve(target.virtual_atoms.size() + source.virtual_atoms.size());
    for (const auto& vatom : source.virtual_atoms) {
        target.virtual_atoms.push_back({vatom.virtual_atom + atom_offset, vatom.atom0 + atom_offset,
                                        vatom.atom1 + atom_offset, vatom.atom2 + atom_offset,
                                        vatom.k1, vatom.k2});
    }
    target.harmonic_impropers.reserve(target.harmonic_impropers.size() + source.harmonic_impropers.size());
    for (const auto& improper : source.harmonic_impropers) {
        target.harmonic_impropers.push_back({improper.atom0 + atom_offset, improper.atom1 + atom_offset,
                                             improper.atom2 + atom_offset, improper.atom3 + atom_offset,
                                             improper.k, improper.phi0});
    }
    const std::uint32_t cmap_type_offset = static_cast<std::uint32_t>(target.cmap_types.size());
    target.cmap_types.reserve(target.cmap_types.size() + source.cmap_types.size());
    for (const auto& type : source.cmap_types) {
        target.cmap_types.push_back(type);
    }
    target.cmaps.reserve(target.cmaps.size() + source.cmaps.size());
    for (const auto& cmap : source.cmaps) {
        target.cmaps.push_back({cmap.atom0 + atom_offset, cmap.atom1 + atom_offset, cmap.atom2 + atom_offset,
                                cmap.atom3 + atom_offset, cmap.atom4 + atom_offset, cmap.type + cmap_type_offset});
    }
    target.nb14_extras.reserve(target.nb14_extras.size() + source.nb14_extras.size());
    for (const auto& nb14 : source.nb14_extras) {
        target.nb14_extras.push_back({nb14.atom1 + atom_offset, nb14.atom2 + atom_offset,
                                      nb14.a, nb14.b, nb14.kee});
    }
    target.urey_bradleys.reserve(target.urey_bradleys.size() + source.urey_bradleys.size());
    for (const auto& angle : source.urey_bradleys) {
        target.urey_bradleys.push_back({angle.atom0 + atom_offset, angle.atom1 + atom_offset,
                                        angle.atom2 + atom_offset, angle.k, angle.b, angle.k_ub, angle.r13});
    }
    target.ryckaert_bellemans.reserve(target.ryckaert_bellemans.size() + source.ryckaert_bellemans.size());
    for (const auto& dihedral : source.ryckaert_bellemans) {
        target.ryckaert_bellemans.push_back({dihedral.atom0 + atom_offset, dihedral.atom1 + atom_offset,
                                             dihedral.atom2 + atom_offset, dihedral.atom3 + atom_offset,
                                             dihedral.c0, dihedral.c1, dihedral.c2,
                                             dihedral.c3, dihedral.c4, dihedral.c5});
    }
    target.soft_bonds.reserve(target.soft_bonds.size() + source.soft_bonds.size());
    for (const auto& bond : source.soft_bonds) {
        target.soft_bonds.push_back({bond.atom1 + atom_offset, bond.atom2 + atom_offset,
                                     bond.k, bond.b, bond.from_a_or_b});
    }
    target.listed_force_definitions.reserve(target.listed_force_definitions.size() +
                                            source.listed_force_definitions.size());
    for (const auto& definition : source.listed_force_definitions) {
        target.listed_force_definitions.push_back(definition);
    }
    for (const auto& [name, parameter] : source.sw_parameters) {
        target.sw_parameters[name] = parameter;
    }
    for (const auto& [name, parameter] : source.edip_parameters) {
        target.edip_parameters[name] = parameter;
    }
}

void remap_internal_structures(const Molecule& source, Molecule& target, const std::vector<AtomId>& old_to_new_atom) {
    target.explicit_bonds.reserve(target.explicit_bonds.size() + source.explicit_bonds.size());
    for (const auto& bond : source.explicit_bonds) {
        const AtomId atom1 = remap_atom_id(old_to_new_atom, bond.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, bond.atom2);
        if (atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        target.explicit_bonds.push_back({atom1, atom2});
    }
    target.residue_links.reserve(target.residue_links.size() + source.residue_links.size());
    for (const auto& link : source.residue_links) {
        const AtomId atom1 = remap_atom_id(old_to_new_atom, link.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, link.atom2);
        if (atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        target.residue_links.push_back({atom1, atom2});
    }
    target.virtual_atoms.reserve(target.virtual_atoms.size() + source.virtual_atoms.size());
    for (const auto& vatom : source.virtual_atoms) {
        const AtomId virtual_atom = remap_atom_id(old_to_new_atom, vatom.virtual_atom);
        const AtomId atom0 = remap_atom_id(old_to_new_atom, vatom.atom0);
        const AtomId atom1 = remap_atom_id(old_to_new_atom, vatom.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, vatom.atom2);
        if (virtual_atom == invalid_atom_id || atom0 == invalid_atom_id ||
            atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        target.virtual_atoms.push_back({virtual_atom, atom0, atom1, atom2, vatom.k1, vatom.k2});
    }
    target.harmonic_impropers.reserve(target.harmonic_impropers.size() + source.harmonic_impropers.size());
    for (const auto& improper : source.harmonic_impropers) {
        const AtomId atom0 = remap_atom_id(old_to_new_atom, improper.atom0);
        const AtomId atom1 = remap_atom_id(old_to_new_atom, improper.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, improper.atom2);
        const AtomId atom3 = remap_atom_id(old_to_new_atom, improper.atom3);
        if (atom0 == invalid_atom_id || atom1 == invalid_atom_id ||
            atom2 == invalid_atom_id || atom3 == invalid_atom_id) {
            continue;
        }
        target.harmonic_impropers.push_back({atom0, atom1, atom2, atom3, improper.k, improper.phi0});
    }
    target.cmap_types = source.cmap_types;
    target.cmaps.reserve(target.cmaps.size() + source.cmaps.size());
    for (const auto& cmap : source.cmaps) {
        const AtomId atom0 = remap_atom_id(old_to_new_atom, cmap.atom0);
        const AtomId atom1 = remap_atom_id(old_to_new_atom, cmap.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, cmap.atom2);
        const AtomId atom3 = remap_atom_id(old_to_new_atom, cmap.atom3);
        const AtomId atom4 = remap_atom_id(old_to_new_atom, cmap.atom4);
        if (atom0 == invalid_atom_id || atom1 == invalid_atom_id || atom2 == invalid_atom_id ||
            atom3 == invalid_atom_id || atom4 == invalid_atom_id) {
            continue;
        }
        target.cmaps.push_back({atom0, atom1, atom2, atom3, atom4, cmap.type});
    }
    target.nb14_extras.reserve(target.nb14_extras.size() + source.nb14_extras.size());
    for (const auto& nb14 : source.nb14_extras) {
        const AtomId atom1 = remap_atom_id(old_to_new_atom, nb14.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, nb14.atom2);
        if (atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        target.nb14_extras.push_back({atom1, atom2, nb14.a, nb14.b, nb14.kee});
    }
    target.urey_bradleys.reserve(target.urey_bradleys.size() + source.urey_bradleys.size());
    for (const auto& angle : source.urey_bradleys) {
        const AtomId atom0 = remap_atom_id(old_to_new_atom, angle.atom0);
        const AtomId atom1 = remap_atom_id(old_to_new_atom, angle.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, angle.atom2);
        if (atom0 == invalid_atom_id || atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        target.urey_bradleys.push_back({atom0, atom1, atom2, angle.k, angle.b, angle.k_ub, angle.r13});
    }
    target.ryckaert_bellemans.reserve(target.ryckaert_bellemans.size() + source.ryckaert_bellemans.size());
    for (const auto& dihedral : source.ryckaert_bellemans) {
        const AtomId atom0 = remap_atom_id(old_to_new_atom, dihedral.atom0);
        const AtomId atom1 = remap_atom_id(old_to_new_atom, dihedral.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, dihedral.atom2);
        const AtomId atom3 = remap_atom_id(old_to_new_atom, dihedral.atom3);
        if (atom0 == invalid_atom_id || atom1 == invalid_atom_id ||
            atom2 == invalid_atom_id || atom3 == invalid_atom_id) {
            continue;
        }
        target.ryckaert_bellemans.push_back({atom0, atom1, atom2, atom3, dihedral.c0, dihedral.c1,
                                             dihedral.c2, dihedral.c3, dihedral.c4, dihedral.c5});
    }
    target.soft_bonds.reserve(target.soft_bonds.size() + source.soft_bonds.size());
    for (const auto& bond : source.soft_bonds) {
        const AtomId atom1 = remap_atom_id(old_to_new_atom, bond.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, bond.atom2);
        if (atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        target.soft_bonds.push_back({atom1, atom2, bond.k, bond.b, bond.from_a_or_b});
    }
}

constexpr std::array<const char*, 112> kElements{
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc",
    "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",
    "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
    "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg",
};

constexpr std::array<double, 112> kElementMasses{
    0.00000, 1.00794, 4.00260, 6.941, 9.012182, 10.811,
    12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
    22.989770, 24.3050, 26.981538, 28.0855, 30.973761,
    32.065, 35.453, 39.948, 39.0983, 40.078, 44.955910,
    47.867, 50.9415, 51.9961, 54.938049, 55.845, 58.9332,
    58.6934, 63.546, 65.409, 69.723, 72.64, 74.92160,
    78.96, 79.904, 83.798, 85.4678, 87.62, 88.90585,
    91.224, 92.90638, 95.94, 98.0, 101.07, 102.90550,
    106.42, 107.8682, 112.411, 114.818, 118.710, 121.760,
    127.60, 126.90447, 131.293, 132.90545, 137.327,
    138.9055, 140.116, 140.90765, 144.24, 145.0, 150.36,
    151.964, 157.25, 158.92534, 162.500, 164.93032,
    167.259, 168.93421, 173.04, 174.967, 178.49, 180.9479,
    183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,
    200.59, 204.3833, 207.2, 208.98038, 209.0, 210.0, 222.0,
    223.0, 226.0, 227.0, 232.0381, 231.03588, 238.02891,
    237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0,
    258.0, 259.0, 262.0, 261.0, 262.0, 266.0, 264.0, 269.0,
    268.0, 271.0, 272.0,
};

std::string fallback_element_from_name(const std::string& atom_name) {
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

std::string canonicalize_explicit_element(const std::string& explicit_element) {
    const auto trimmed = trim_copy(explicit_element);
    if (trimmed.empty()) {
        return "";
    }
    std::string letters;
    for (const char c : trimmed) {
        if (std::isalpha(static_cast<unsigned char>(c))) {
            letters.push_back(c);
        }
    }
    if (letters.empty()) {
        return "";
    }
    std::string element;
    element.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(letters[0]))));
    if (letters.size() > 1) {
        element.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(letters[1]))));
    }
    return element;
}

std::optional<double> registered_mass_for_type(const std::string& atom_type) {
    if (atom_type.empty()) {
        return std::nullopt;
    }
    if (const auto amber_mass = find_amber_atom_type_mass(atom_type)) {
        return amber_mass;
    }
    if (const auto external_mass = find_external_atom_type_mass(atom_type)) {
        return external_mass;
    }
    return std::nullopt;
}

double resolve_mass_for_atom(const std::string& atom_name, const std::string& atom_type,
                             const std::string& explicit_element, double explicit_mass) {
    if (explicit_mass > 0.0) {
        return explicit_mass;
    }
    if (const auto registered_mass = registered_mass_for_type(atom_type)) {
        return *registered_mass;
    }
    const auto explicit_canonical = canonicalize_explicit_element(explicit_element);
    if (!explicit_canonical.empty()) {
        const auto default_mass = default_mass_for_element(explicit_canonical);
        if (default_mass > 0.0) {
            return default_mass;
        }
    }
    const auto fallback_element = fallback_element_from_name(atom_name);
    return default_mass_for_element(fallback_element);
}

std::string resolve_element_for_atom(const std::string& atom_name, const std::string& explicit_element,
                                     double resolved_mass) {
    const auto explicit_canonical = canonicalize_explicit_element(explicit_element);
    if (!explicit_canonical.empty()) {
        return explicit_canonical;
    }
    if (resolved_mass > 0.0) {
        return guess_element_from_mass(resolved_mass);
    }
    return fallback_element_from_name(atom_name);
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

std::string guess_element_from_mass(double mass) {
    std::size_t index = 0;
    if (mass > 0.0 && mass < 3.8) {
        index = 1;
    } else if (mass > 207.85 && mass < 208.99) {
        index = 83;
    } else if (mass > 56.50 && mass < 58.8133) {
        index = 27;
    } else {
        for (std::size_t j = 0; j < kElementMasses.size(); ++j) {
            if (std::abs(mass - kElementMasses[j]) < 0.65) {
                index = j;
                break;
            }
        }
    }
    return kElements[index];
}

std::string guess_element(const std::string& atom_name, const std::string& explicit_element) {
    const auto explicit_canonical = canonicalize_explicit_element(explicit_element);
    if (!explicit_canonical.empty()) {
        return explicit_canonical;
    }
    return fallback_element_from_name(atom_name);
}

ResidueType::ResidueType(std::string name) : name_(std::move(name)) {}

const std::string& ResidueType::name() const noexcept { return name_; }
std::uint64_t ResidueType::version() const noexcept { return version_; }
std::size_t ResidueType::atom_count() const noexcept { return atoms_.size(); }
std::size_t ResidueType::bond_count() const noexcept { return bonds_.size(); }
const std::vector<ResidueTypeAtom>& ResidueType::atoms() const noexcept { return atoms_; }
const std::vector<ResidueTypeBond>& ResidueType::bonds() const noexcept { return bonds_; }
const std::string& ResidueType::head() const noexcept { return head_; }
const std::string& ResidueType::tail() const noexcept { return tail_; }
const std::string& ResidueType::head_next() const noexcept { return head_next_; }
const std::string& ResidueType::tail_next() const noexcept { return tail_next_; }
double ResidueType::head_length() const noexcept { return head_length_; }
double ResidueType::tail_length() const noexcept { return tail_length_; }
const std::unordered_map<std::string, std::string>& ResidueType::connect_atoms() const noexcept {
    return connect_atoms_;
}

void ResidueType::add_atom(const std::string& name, const std::string& type, double x, double y, double z,
                           double charge, double mass) {
    if (atom_name_to_index_.count(name) != 0) {
        throw std::invalid_argument("duplicate atom name in ResidueType: " + name);
    }
    ResidueTypeAtom atom;
    atom.name = name;
    atom.type = type;
    atom.mass = resolve_mass_for_atom(name, type, "", mass);
    atom.element = resolve_element_for_atom(name, "", atom.mass);
    atom.x = x;
    atom.y = y;
    atom.z = z;
    atom.charge = charge;
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

void ResidueType::set_head(const std::string& atom, double length, const std::string& next) {
    head_ = atom;
    head_length_ = length;
    head_next_ = next;
    ++version_;
}

void ResidueType::set_tail(const std::string& atom, double length, const std::string& next) {
    tail_ = atom;
    tail_length_ = length;
    tail_next_ = next;
    ++version_;
}

void ResidueType::set_connect_atom(const std::string& key, const std::string& atom) {
    connect_atoms_[key] = atom;
    ++version_;
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
    residue.original_name = type.name();
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

void Molecule::add_molecule(const Molecule& other) {
    add_molecule_linked(other, false);
}

void Molecule::add_molecule_linked(const Molecule& other, bool link) {
    if (!other.validate()) {
        throw std::invalid_argument("source molecule is invalid");
    }
    std::optional<AtomId> link_atom1;
    std::optional<AtomId> link_atom2;
    if (link && !residues.empty() && !other.residues.empty()) {
        const auto& left = residues.back();
        const auto& right = other.residues.front();
        if (has_template(left.name) && has_template(right.name)) {
            const auto& left_type = get_residue_template(left.name);
            const auto& right_type = get_residue_template(right.name);
            if (!left_type.tail().empty() && !right_type.head().empty()) {
                for (std::uint32_t local = 0; local < left.atom_count; ++local) {
                    const AtomId atom_id = left.atom_begin + local;
                    if (atoms[atom_id].name == left_type.tail()) {
                        link_atom1 = atom_id;
                        break;
                    }
                }
                for (std::uint32_t local = 0; local < right.atom_count; ++local) {
                    const AtomId atom_id = right.atom_begin + local;
                    if (other.atoms[atom_id].name == right_type.head()) {
                        link_atom2 = atom_id;
                        break;
                    }
                }
            }
        }
    }
    const AtomId atom_offset = static_cast<AtomId>(atoms.size());
    const ResidueId residue_offset = static_cast<ResidueId>(residues.size());

    atoms.reserve(atoms.size() + other.atoms.size());
    residues.reserve(residues.size() + other.residues.size());
    explicit_bonds.reserve(explicit_bonds.size() + other.explicit_bonds.size());
    residue_links.reserve(residue_links.size() + other.residue_links.size());
    virtual_atoms.reserve(virtual_atoms.size() + other.virtual_atoms.size());
    harmonic_impropers.reserve(harmonic_impropers.size() + other.harmonic_impropers.size());
    cmap_types.reserve(cmap_types.size() + other.cmap_types.size());
    cmaps.reserve(cmaps.size() + other.cmaps.size());
    nb14_extras.reserve(nb14_extras.size() + other.nb14_extras.size());
    urey_bradleys.reserve(urey_bradleys.size() + other.urey_bradleys.size());
    ryckaert_bellemans.reserve(ryckaert_bellemans.size() + other.ryckaert_bellemans.size());
    soft_bonds.reserve(soft_bonds.size() + other.soft_bonds.size());
    listed_force_definitions.reserve(listed_force_definitions.size() + other.listed_force_definitions.size());

    for (const auto& residue : other.residues) {
        Residue copied = residue;
        copied.atom_begin += atom_offset;
        residues.push_back(std::move(copied));
    }

    for (const auto& source_atom : other.atoms) {
        Atom atom = source_atom;
        atom.residue += residue_offset;
        atoms.push_back(std::move(atom));
    }

    for (const auto& bond : other.explicit_bonds) {
        explicit_bonds.push_back({bond.atom1 + atom_offset, bond.atom2 + atom_offset});
    }
    for (const auto& link : other.residue_links) {
        residue_links.push_back({link.atom1 + atom_offset, link.atom2 + atom_offset});
    }
    for (const auto& vatom : other.virtual_atoms) {
        virtual_atoms.push_back({vatom.virtual_atom + atom_offset, vatom.atom0 + atom_offset,
                                 vatom.atom1 + atom_offset, vatom.atom2 + atom_offset,
                                 vatom.k1, vatom.k2});
    }
    for (const auto& improper : other.harmonic_impropers) {
        harmonic_impropers.push_back({improper.atom0 + atom_offset, improper.atom1 + atom_offset,
                                      improper.atom2 + atom_offset, improper.atom3 + atom_offset,
                                      improper.k, improper.phi0});
    }
    const std::uint32_t cmap_type_offset = static_cast<std::uint32_t>(cmap_types.size());
    for (const auto& type : other.cmap_types) {
        cmap_types.push_back(type);
    }
    for (const auto& cmap : other.cmaps) {
        cmaps.push_back({cmap.atom0 + atom_offset, cmap.atom1 + atom_offset, cmap.atom2 + atom_offset,
                         cmap.atom3 + atom_offset, cmap.atom4 + atom_offset, cmap.type + cmap_type_offset});
    }
    for (const auto& nb14 : other.nb14_extras) {
        nb14_extras.push_back({nb14.atom1 + atom_offset, nb14.atom2 + atom_offset, nb14.a, nb14.b, nb14.kee});
    }
    for (const auto& angle : other.urey_bradleys) {
        urey_bradleys.push_back({angle.atom0 + atom_offset, angle.atom1 + atom_offset, angle.atom2 + atom_offset,
                                 angle.k, angle.b, angle.k_ub, angle.r13});
    }
    for (const auto& dihedral : other.ryckaert_bellemans) {
        ryckaert_bellemans.push_back({dihedral.atom0 + atom_offset, dihedral.atom1 + atom_offset,
                                      dihedral.atom2 + atom_offset, dihedral.atom3 + atom_offset,
                                      dihedral.c0, dihedral.c1, dihedral.c2,
                                      dihedral.c3, dihedral.c4, dihedral.c5});
    }
    for (const auto& bond : other.soft_bonds) {
        soft_bonds.push_back({bond.atom1 + atom_offset, bond.atom2 + atom_offset, bond.k, bond.b,
                              bond.from_a_or_b});
    }
    for (const auto& definition : other.listed_force_definitions) {
        listed_force_definitions.push_back(definition);
    }
    for (const auto& [name, parameter] : other.sw_parameters) {
        sw_parameters[name] = parameter;
    }
    for (const auto& [name, parameter] : other.edip_parameters) {
        edip_parameters[name] = parameter;
    }

    if (!validate()) {
        throw std::runtime_error("internal error: invalid molecule after merge");
    }
    if (link_atom1 && link_atom2) {
        add_residue_link(*link_atom1, *link_atom2 + atom_offset);
    }
}

void Molecule::add_residue_link(AtomId atom1, AtomId atom2) {
    ensure_atom_id(*this, atom1);
    ensure_atom_id(*this, atom2);
    if (atom1 == atom2) {
        return;
    }
    if (atoms[atom1].residue == atoms[atom2].residue) {
        return;
    }
    const auto lo = std::min(atom1, atom2);
    const auto hi = std::max(atom1, atom2);
    for (const auto& link : residue_links) {
        if (std::min(link.atom1, link.atom2) == lo && std::max(link.atom1, link.atom2) == hi) {
            return;
        }
    }
    residue_links.push_back({lo, hi});
}

void Molecule::add_virtual_atom2(AtomId virtual_atom, AtomId atom0, AtomId atom1, AtomId atom2,
                                 double k1, double k2) {
    ensure_atom_id(*this, virtual_atom);
    ensure_atom_id(*this, atom0);
    ensure_atom_id(*this, atom1);
    ensure_atom_id(*this, atom2);
    virtual_atoms.push_back({virtual_atom, atom0, atom1, atom2, k1, k2});
}

void Molecule::add_improper_dihedral(AtomId atom0, AtomId atom1, AtomId atom2, AtomId atom3,
                                     double k, double phi0) {
    ensure_atom_id(*this, atom0);
    ensure_atom_id(*this, atom1);
    ensure_atom_id(*this, atom2);
    ensure_atom_id(*this, atom3);
    harmonic_impropers.push_back({atom0, atom1, atom2, atom3, k, phi0});
}

std::uint32_t Molecule::add_cmap_type(std::uint32_t resolution, const std::vector<double>& parameters) {
    if (resolution == 0) {
        throw std::invalid_argument("cmap resolution should be positive");
    }
    const auto expected = static_cast<std::size_t>(resolution) * resolution;
    if (parameters.size() != expected) {
        throw std::invalid_argument("cmap parameter count should equal resolution * resolution");
    }
    const auto index = static_cast<std::uint32_t>(cmap_types.size());
    cmap_types.push_back({resolution, parameters});
    return index;
}

void Molecule::add_cmap(AtomId atom0, AtomId atom1, AtomId atom2, AtomId atom3, AtomId atom4,
                        std::uint32_t type) {
    ensure_atom_id(*this, atom0);
    ensure_atom_id(*this, atom1);
    ensure_atom_id(*this, atom2);
    ensure_atom_id(*this, atom3);
    ensure_atom_id(*this, atom4);
    if (type >= cmap_types.size()) {
        throw std::out_of_range("cmap type out of range");
    }
    cmaps.push_back({atom0, atom1, atom2, atom3, atom4, type});
}

void Molecule::add_nb14_extra(AtomId atom1, AtomId atom2, double a, double b, double kee) {
    ensure_atom_id(*this, atom1);
    ensure_atom_id(*this, atom2);
    if (atom1 == atom2) {
        throw std::invalid_argument("nb14_extra atoms should be different");
    }
    nb14_extras.push_back({atom1, atom2, a, b, kee});
}

void Molecule::add_urey_bradley(AtomId atom0, AtomId atom1, AtomId atom2,
                                double k, double b, double k_ub, double r13) {
    ensure_atom_id(*this, atom0);
    ensure_atom_id(*this, atom1);
    ensure_atom_id(*this, atom2);
    urey_bradleys.push_back({atom0, atom1, atom2, k, b, k_ub, r13});
}

void Molecule::add_ryckaert_bellemans(AtomId atom0, AtomId atom1, AtomId atom2, AtomId atom3,
                                      double c0, double c1, double c2, double c3, double c4, double c5) {
    ensure_atom_id(*this, atom0);
    ensure_atom_id(*this, atom1);
    ensure_atom_id(*this, atom2);
    ensure_atom_id(*this, atom3);
    ryckaert_bellemans.push_back({atom0, atom1, atom2, atom3, c0, c1, c2, c3, c4, c5});
}

void Molecule::add_bond_soft(AtomId atom1, AtomId atom2, double k, double b, int from_a_or_b) {
    ensure_atom_id(*this, atom1);
    ensure_atom_id(*this, atom2);
    if (atom1 == atom2) {
        throw std::invalid_argument("bond_soft atoms should be different");
    }
    soft_bonds.push_back({atom1, atom2, k, b, from_a_or_b});
}

void Molecule::add_listed_force_definition(const std::string& definition) {
    if (!definition.empty()) {
        listed_force_definitions.push_back(definition);
    }
}

namespace {

std::pair<double, double> bondi_radius_and_scaler(const std::string& element) {
    if (element == "H") return {1.2, 0.85};
    if (element == "C") return {1.7, 0.72};
    if (element == "N") return {1.55, 0.79};
    if (element == "O") return {1.52, 0.85};
    if (element == "F") return {1.47, 0.88};
    if (element == "P") return {1.8, 0.86};
    if (element == "S") return {1.8, 0.96};
    if (element == "Cl") return {1.75, 0.8};
    if (element == "Br") return {1.85, 0.8};
    if (element == "I") return {1.98, 0.8};
    return {1.5, 0.8};
}

std::pair<double, double> modified_bondi_radius_and_scaler(const std::string& element) {
    if (element == "H") return {1.2, 0.85};
    if (element == "C") return {1.7, 0.72};
    if (element == "N") return {1.55, 0.79};
    if (element == "O") return {1.5, 0.85};
    if (element == "F") return {1.5, 0.88};
    if (element == "Si") return {2.1, 0.8};
    if (element == "P") return {1.85, 0.86};
    if (element == "S") return {1.8, 0.96};
    if (element == "Cl") return {1.7, 0.8};
    if (element == "Br") return {1.85, 0.8};
    if (element == "I") return {1.98, 0.8};
    return {1.5, 0.8};
}

std::string gb_element_from_mass(const Atom& atom) {
    return atom.mass > 0.0 ? guess_element_from_mass(atom.mass) : atom.element;
}

std::string gb_element_from_mass_value(double mass, const std::string& fallback_element = "") {
    return mass > 0.0 ? guess_element_from_mass(mass) : fallback_element;
}

}  // namespace

void Molecule::set_gb_radius(const std::string& radius_set) {
    for (AtomId atom_id = 0; atom_id < atoms.size(); ++atom_id) {
        auto& atom = atoms[atom_id];
        const auto element = gb_element_from_mass(atom);
        if (radius_set == "bondi_radii") {
            const auto [radius, scaler] = bondi_radius_and_scaler(element);
            atom.gb_radius = radius;
            atom.gb_scaler = scaler;
        } else if (radius_set == "modified_bondi_radii") {
            auto [radius, scaler] = modified_bondi_radius_and_scaler(element);
            if (element == "H") {
                AtomId bonded_atom_id = invalid_atom_id;
                for (const auto& bond : explicit_bonds) {
                    if (bond.atom1 == atom_id) {
                        bonded_atom_id = bond.atom2;
                        break;
                    }
                    if (bond.atom2 == atom_id) {
                        bonded_atom_id = bond.atom1;
                        break;
                    }
                }
                const auto residue_id = atom.residue;
                if (bonded_atom_id != invalid_atom_id && bonded_atom_id < atoms.size()) {
                    const auto bonded_element = gb_element_from_mass(atoms[bonded_atom_id]);
                    if (bonded_element == "C" || bonded_element == "N") {
                        radius = 1.3;
                    } else if (bonded_element == "S" || bonded_element == "O" || bonded_element == "H") {
                        radius = 0.8;
                    }
                } else if (residue_id < residues.size()) {
                    const auto& residue = residues[residue_id];
                    if (has_template(residue.name)) {
                        const auto& residue_type = get_residue_template(residue.name);
                        auto template_atom_index = invalid_atom_id;
                        try {
                            template_atom_index = residue_type.atom_index(atom.name);
                        } catch (const std::exception&) {
                            template_atom_index = invalid_atom_id;
                        }
                        if (template_atom_index != invalid_atom_id) {
                            for (const auto& bond : residue_type.bonds()) {
                                AtomId bonded_template = invalid_atom_id;
                                if (bond.atom1 == template_atom_index) {
                                    bonded_template = bond.atom2;
                                } else if (bond.atom2 == template_atom_index) {
                                    bonded_template = bond.atom1;
                                }
                                if (bonded_template == invalid_atom_id) {
                                    continue;
                                }
                                const auto& bonded_atom = residue_type.atoms()[bonded_template];
                                const auto bonded_element =
                                    !bonded_atom.element.empty() ? bonded_atom.element :
                                                                  gb_element_from_mass_value(bonded_atom.mass);
                                if (bonded_element == "C" || bonded_element == "N") {
                                    radius = 1.3;
                                } else if (bonded_element == "S" || bonded_element == "O" || bonded_element == "H") {
                                    radius = 0.8;
                                }
                                break;
                            }
                        }
                    }
                }
            }
            atom.gb_radius = radius;
            atom.gb_scaler = scaler;
        } else {
            throw std::invalid_argument("unknown GB radius set: " + radius_set);
        }
    }
    box_length = {999.0, 999.0, 999.0};
    has_box = true;
    has_gb_parameters = true;
}

void Molecule::enable_min_bonded_parameters(bool enabled) noexcept {
    write_min_bonded_parameters = enabled;
}

void Molecule::enable_subsys_division(bool enabled) noexcept {
    write_subsys_division = enabled;
}

void Molecule::enable_lj_soft_core(bool enabled) noexcept {
    write_lj_soft_core = enabled;
    if (enabled) {
        write_subsys_division = true;
    }
}

void Molecule::add_sw_type(const std::string& name, double a_big, double b_big, double epsilon, double p, double q,
                           double a, double gamma, double sigma, double lambda, double b) {
    if (name.empty()) {
        throw std::invalid_argument("SW type name should not be empty");
    }
    sw_parameters[name] = {a_big, b_big, epsilon, p, q, a, gamma, sigma, lambda, b};
}

void Molecule::add_edip_type(const std::string& name, double a_big, double b_big, double a, double c, double alpha,
                             double beta, double eta, double gamma, double lambda, double mu, double rho,
                             double sigma, double q0, double u1, double u2, double u3, double u4) {
    if (name.empty()) {
        throw std::invalid_argument("EDIP type name should not be empty");
    }
    edip_parameters[name] = {a_big, b_big, a, c, alpha, beta, eta, gamma, lambda, mu, rho, sigma,
                             q0, u1, u2, u3, u4};
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

void Molecule::replace_residues(const std::unordered_map<ResidueId, Molecule>& replacements,
                                const std::vector<double>& residue_sort_keys, bool sort) {
    for (const auto& [residue_id, replacement] : replacements) {
        ensure_residue_id(*this, residue_id);
        if (replacement.residue_count() != 1) {
            throw std::invalid_argument("replacement molecules should contain exactly one residue");
        }
        if (!replacement.validate()) {
            throw std::invalid_argument("replacement molecule is invalid");
        }
    }
    if (!residue_sort_keys.empty() && residue_sort_keys.size() != residues.size()) {
        throw std::invalid_argument("residue_sort_keys should match residue count");
    }

    Molecule rebuilt(name);
    rebuilt.box_length = box_length;
    rebuilt.box_angle = box_angle;
    rebuilt.has_box = has_box;
    rebuilt.has_gb_parameters = has_gb_parameters;
    rebuilt.write_min_bonded_parameters = write_min_bonded_parameters;
    rebuilt.write_subsys_division = write_subsys_division;
    rebuilt.write_lj_soft_core = write_lj_soft_core;
    rebuilt.sw_parameters = sw_parameters;
    rebuilt.edip_parameters = edip_parameters;

    std::vector<AtomId> old_to_new_atom(atoms.size(), invalid_atom_id);
    std::vector<ResidueId> residue_order(residues.size());
    for (ResidueId residue_id = 0; residue_id < residues.size(); ++residue_id) {
        residue_order[residue_id] = residue_id;
    }
    if (sort && !residue_sort_keys.empty()) {
        std::stable_sort(residue_order.begin(), residue_order.end(),
                         [&](ResidueId lhs, ResidueId rhs) {
                             return residue_sort_keys[lhs] < residue_sort_keys[rhs];
                         });
    }

    rebuilt.residues.reserve(residues.size() - replacements.size() + replacements.size());
    for (const ResidueId old_residue_id : residue_order) {
        const auto replacement_it = replacements.find(old_residue_id);
        if (replacement_it == replacements.end()) {
            const auto& source_residue = residues[old_residue_id];
            const ResidueId new_residue_id = static_cast<ResidueId>(rebuilt.residues.size());
            Residue copied_residue = source_residue;
            copied_residue.atom_begin = static_cast<AtomId>(rebuilt.atoms.size());
            rebuilt.residues.push_back(copied_residue);
            for (std::uint32_t local = 0; local < source_residue.atom_count; ++local) {
                const AtomId old_atom_id = source_residue.atom_begin + local;
                Atom atom_copy = atoms[old_atom_id];
                atom_copy.residue = new_residue_id;
                const AtomId new_atom_id = static_cast<AtomId>(rebuilt.atoms.size());
                rebuilt.atoms.push_back(std::move(atom_copy));
                old_to_new_atom[old_atom_id] = new_atom_id;
            }
            continue;
        }

        const auto& source_residue = residues[old_residue_id];
        const auto& replacement = replacement_it->second;
        const auto& replacement_residue = replacement.residues.front();
        const ResidueId new_residue_id = static_cast<ResidueId>(rebuilt.residues.size());
        const AtomId atom_offset = static_cast<AtomId>(rebuilt.atoms.size());
        Residue copied_residue = replacement_residue;
        copied_residue.chain_id = source_residue.chain_id;
        copied_residue.effective_chain_id = source_residue.effective_chain_id;
        copied_residue.segment_id = source_residue.segment_id;
        copied_residue.pdb_resseq = source_residue.pdb_resseq;
        copied_residue.insertion_code = source_residue.insertion_code;
        copied_residue.is_hetero = source_residue.is_hetero;
        copied_residue.atom_begin = atom_offset;
        copied_residue.atom_count = replacement_residue.atom_count;
        rebuilt.residues.push_back(copied_residue);

        double dx = 0.0;
        double dy = 0.0;
        double dz = 0.0;
        if (source_residue.atom_count > 0 && replacement_residue.atom_count > 0) {
            const auto& source_anchor = atoms[source_residue.atom_begin];
            const auto& replacement_anchor = replacement.atoms[replacement_residue.atom_begin];
            dx = source_anchor.x - replacement_anchor.x;
            dy = source_anchor.y - replacement_anchor.y;
            dz = source_anchor.z - replacement_anchor.z;
        }

        for (std::uint32_t local = 0; local < replacement_residue.atom_count; ++local) {
            Atom atom_copy = replacement.atoms[replacement_residue.atom_begin + local];
            atom_copy.residue = new_residue_id;
            atom_copy.x += dx;
            atom_copy.y += dy;
            atom_copy.z += dz;
            rebuilt.atoms.push_back(std::move(atom_copy));
        }
        append_internal_structures(rebuilt, replacement, atom_offset);
    }

    remap_internal_structures(*this, rebuilt, old_to_new_atom);
    *this = std::move(rebuilt);
    if (!validate()) {
        throw std::runtime_error("internal error: invalid molecule after residue replacement");
    }
}

void Molecule::reorder_atoms_by_template(const Molecule& template_molecule) {
    if (residue_count() != template_molecule.residue_count()) {
        throw std::invalid_argument("template molecule should have the same residue count");
    }

    Molecule rebuilt(name);
    rebuilt.box_length = box_length;
    rebuilt.box_angle = box_angle;
    rebuilt.has_box = has_box;
    rebuilt.has_gb_parameters = has_gb_parameters;
    rebuilt.write_min_bonded_parameters = write_min_bonded_parameters;
    rebuilt.write_subsys_division = write_subsys_division;
    rebuilt.write_lj_soft_core = write_lj_soft_core;

    std::vector<AtomId> old_to_new_atom(atoms.size(), invalid_atom_id);
    rebuilt.residues.reserve(residues.size());
    rebuilt.atoms.reserve(atoms.size());

    for (ResidueId residue_id = 0; residue_id < residues.size(); ++residue_id) {
        const auto& source_residue = residues[residue_id];
        const auto& template_residue = template_molecule.residues[residue_id];
        const auto& source_type = !source_residue.type_name.empty() ? source_residue.type_name : source_residue.name;
        const auto& template_type = !template_residue.type_name.empty() ? template_residue.type_name : template_residue.name;
        if (source_type != template_type) {
            throw std::invalid_argument("residue types should match when sorting atoms by template");
        }
        if (source_residue.atom_count != template_residue.atom_count) {
            throw std::invalid_argument("residue atom counts should match when sorting atoms by template");
        }

        std::unordered_map<std::string, AtomId> source_name_to_atom;
        source_name_to_atom.reserve(source_residue.atom_count);
        for (std::uint32_t local = 0; local < source_residue.atom_count; ++local) {
            const AtomId atom_id = source_residue.atom_begin + local;
            const auto inserted = source_name_to_atom.emplace(atoms[atom_id].name, atom_id);
            if (!inserted.second) {
                throw std::invalid_argument("duplicate atom names are not supported in sort_atoms_by");
            }
        }

        Residue copied_residue = source_residue;
        copied_residue.atom_begin = static_cast<AtomId>(rebuilt.atoms.size());
        rebuilt.residues.push_back(copied_residue);
        for (std::uint32_t local = 0; local < template_residue.atom_count; ++local) {
            const auto& template_atom = template_molecule.atoms[template_residue.atom_begin + local];
            const auto found = source_name_to_atom.find(template_atom.name);
            if (found == source_name_to_atom.end()) {
                throw std::invalid_argument("template atom name not found in source residue: " + template_atom.name);
            }
            const AtomId old_atom_id = found->second;
            Atom atom_copy = atoms[old_atom_id];
            atom_copy.residue = residue_id;
            const AtomId new_atom_id = static_cast<AtomId>(rebuilt.atoms.size());
            rebuilt.atoms.push_back(std::move(atom_copy));
            old_to_new_atom[old_atom_id] = new_atom_id;
        }
    }

    remap_internal_structures(*this, rebuilt, old_to_new_atom);
    rebuilt.listed_force_definitions = listed_force_definitions;
    rebuilt.sw_parameters = sw_parameters;
    rebuilt.edip_parameters = edip_parameters;
    *this = std::move(rebuilt);
    if (!validate()) {
        throw std::runtime_error("internal error: invalid molecule after atom reordering");
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
    for (const auto& vatom : virtual_atoms) {
        if (vatom.virtual_atom >= atoms.size() || vatom.atom0 >= atoms.size() ||
            vatom.atom1 >= atoms.size() || vatom.atom2 >= atoms.size()) {
            return false;
        }
    }
    for (const auto& improper : harmonic_impropers) {
        if (improper.atom0 >= atoms.size() || improper.atom1 >= atoms.size() ||
            improper.atom2 >= atoms.size() || improper.atom3 >= atoms.size()) {
            return false;
        }
    }
    for (const auto& type : cmap_types) {
        if (type.resolution == 0 ||
            type.parameters.size() != static_cast<std::size_t>(type.resolution) * type.resolution) {
            return false;
        }
    }
    for (const auto& cmap : cmaps) {
        if (cmap.atom0 >= atoms.size() || cmap.atom1 >= atoms.size() || cmap.atom2 >= atoms.size() ||
            cmap.atom3 >= atoms.size() || cmap.atom4 >= atoms.size() || cmap.type >= cmap_types.size()) {
            return false;
        }
    }
    for (const auto& nb14 : nb14_extras) {
        if (nb14.atom1 >= atoms.size() || nb14.atom2 >= atoms.size() || nb14.atom1 == nb14.atom2) {
            return false;
        }
    }
    for (const auto& angle : urey_bradleys) {
        if (angle.atom0 >= atoms.size() || angle.atom1 >= atoms.size() || angle.atom2 >= atoms.size()) {
            return false;
        }
    }
    for (const auto& dihedral : ryckaert_bellemans) {
        if (dihedral.atom0 >= atoms.size() || dihedral.atom1 >= atoms.size() ||
            dihedral.atom2 >= atoms.size() || dihedral.atom3 >= atoms.size()) {
            return false;
        }
    }
    for (const auto& bond : soft_bonds) {
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

}  // namespace xpongecpp
