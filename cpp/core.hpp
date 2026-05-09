#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace xpongecpp {

using AtomId = std::uint32_t;
using ResidueId = std::uint32_t;

struct Atom {
    std::string name;
    std::string type;
    std::string element;
    ResidueId residue{0};
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double charge{0.0};
    double mass{0.0};
};

struct Residue {
    std::string name;
    std::string type_name;
    AtomId atom_begin{0};
    std::uint32_t atom_count{0};
};

struct ResidueLink {
    AtomId atom1{0};
    AtomId atom2{0};
};

struct Bond {
    AtomId atom1{0};
    AtomId atom2{0};
    double k{0.0};
    double length{0.0};
};

struct Angle {
    AtomId atom1{0};
    AtomId atom2{0};
    AtomId atom3{0};
    double k{0.0};
    double theta{0.0};
};

struct Dihedral {
    AtomId atom1{0};
    AtomId atom2{0};
    AtomId atom3{0};
    AtomId atom4{0};
    int periodicity{1};
    double k{0.0};
    double phase{0.0};
};

struct NB14 {
    AtomId atom1{0};
    AtomId atom2{0};
    double k_lj{0.5};
    double k_ee{0.833333};
};

struct Topology {
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;
    std::vector<std::vector<AtomId>> exclusions;
    std::vector<NB14> nb14s;
};

struct ResidueTypeAtom {
    std::string name;
    std::string type;
    std::string element;
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double charge{0.0};
    double mass{0.0};
};

struct ResidueTypeBond {
    std::uint32_t atom1{0};
    std::uint32_t atom2{0};
};

class ResidueType {
public:
    explicit ResidueType(std::string name);

    const std::string& name() const noexcept;
    std::uint64_t version() const noexcept;
    std::size_t atom_count() const noexcept;
    std::size_t bond_count() const noexcept;
    const std::vector<ResidueTypeAtom>& atoms() const noexcept;
    const std::vector<ResidueTypeBond>& bonds() const noexcept;

    void add_atom(const std::string& name, const std::string& type, double x, double y, double z,
                  double charge, double mass);
    void add_connectivity(const std::string& atom1, const std::string& atom2);
    std::uint32_t atom_index(const std::string& name) const;

private:
    std::string name_;
    std::uint64_t version_{0};
    std::vector<ResidueTypeAtom> atoms_;
    std::vector<ResidueTypeBond> bonds_;
    std::unordered_map<std::string, std::uint32_t> atom_name_to_index_;
};

class Molecule {
public:
    explicit Molecule(std::string name = "MOL");

    std::string name;
    std::vector<Atom> atoms;
    std::vector<Residue> residues;
    std::vector<ResidueLink> residue_links;
    std::array<double, 3> box_length{0.0, 0.0, 0.0};
    std::array<double, 3> box_angle{90.0, 90.0, 90.0};

    std::size_t atom_count() const noexcept;
    std::size_t residue_count() const noexcept;
    const Atom& atom(AtomId id) const;
    Atom& atom(AtomId id);
    const Residue& residue(ResidueId id) const;
    Residue& residue(ResidueId id);
    void append_residue_from_type(const ResidueType& type, double dx, double dy, double dz);
    void set_box_padding(double padding, bool center);
    bool validate() const;
    std::unordered_map<std::string, std::size_t> residue_counts() const;
};

class Assign {
public:
    explicit Assign(std::string name = "ASN");

    std::string name;
    std::vector<std::string> elements;
    std::vector<std::string> names;
    std::vector<std::array<double, 3>> coordinates;
    std::vector<double> charges;
    std::vector<std::unordered_map<std::uint32_t, int>> bonds;
    std::vector<std::string> atom_types;

    void add_atom(const std::string& element, double x, double y, double z,
                  const std::string& name = "", double charge = 0.0);
    void add_bond(std::uint32_t atom1, std::uint32_t atom2, int order);
    void determine_connectivity(double simple_cutoff);
    void determine_atom_type(const std::string& rule);
    ResidueType to_residuetype(const std::string& name) const;
};

Molecule load_pdb_text(const std::string& text);
Molecule load_mol2_text(const std::string& text);
void add_solvent_box(Molecule& molecule, const Molecule& solvent, double distance, double tolerance,
                     std::int64_t n_solvent);
void add_ions(Molecule& molecule, const std::unordered_map<std::string, std::int64_t>& counts);
std::unordered_map<std::string, std::filesystem::path> save_sponge_input(const Molecule& molecule,
                                                                         const std::string& prefix,
                                                                         const std::filesystem::path& dirname);
void save_pdb(const Molecule& molecule, const std::filesystem::path& filename);

Topology build_topology(const Molecule& molecule);

void register_ff14sb();
void register_tip3p();
void register_residue_templates_from_mol2_text(const std::string& text);
void register_residue_templates_from_mol2_file(const std::filesystem::path& filename);
bool has_template(const std::string& name);
std::size_t template_atom_count(const std::string& name);
Molecule get_template_molecule(const std::string& name);
const ResidueType& get_residue_template(const std::string& name);

double default_mass_for_element(const std::string& element);
std::string guess_element(const std::string& atom_name, const std::string& explicit_element = "");

}  // namespace xpongecpp
