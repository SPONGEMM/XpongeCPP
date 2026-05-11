#include "core.hpp"

#include <algorithm>
#include <cctype>
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

}  // namespace

void Molecule::set_gb_radius(const std::string& radius_set) {
    for (auto& atom : atoms) {
        if (radius_set == "bondi_radii") {
            const auto [radius, scaler] = bondi_radius_and_scaler(atom.element);
            atom.gb_radius = radius;
            atom.gb_scaler = scaler;
        } else if (radius_set == "modified_bondi_radii") {
            const auto [radius, scaler] = modified_bondi_radius_and_scaler(atom.element);
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
