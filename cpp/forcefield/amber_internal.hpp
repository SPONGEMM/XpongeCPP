#pragma once

#include "core.hpp"

#include <array>
#include <shared_mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace xpongecpp {

struct ProperParameter {
    std::array<std::string, 4> types;
    std::vector<DihedralTerm> terms;
    std::size_t order{0};
};

struct ImproperParameter {
    std::array<std::string, 4> types;
    DihedralTerm term;
    std::size_t order{0};
};

struct NB14Parameter {
    std::string atom_type1;
    std::string atom_type4;
    NB14Scale scale;
    std::size_t order{0};
};

struct AmberCMapParameter {
    std::uint32_t resolution{0};
    std::vector<double> parameters;
};

struct BondParameter {
    std::string atom_type1;
    std::string atom_type2;
    BondTerm term;
    std::size_t order{0};
};

struct AngleParameter {
    std::array<std::string, 3> types;
    AngleTerm term;
    std::size_t order{0};
};

std::unordered_map<std::string, ResidueType>& templates();
std::unordered_map<std::string, Molecule>& molecule_templates();
std::shared_mutex& registry_mutex();

std::vector<BondParameter>& bond_parameters();
std::vector<AngleParameter>& angle_parameters();
std::vector<ProperParameter>& proper_parameters();
std::vector<ImproperParameter>& improper_parameters();
std::vector<NB14Parameter>& nb14_parameters();
std::unordered_map<std::string, AmberCMapParameter>& amber_cmap_parameters();
std::unordered_map<std::string, std::string>& lj_type_by_atom_type();
std::unordered_map<std::string, std::pair<double, double>>& lj_parameters();
std::unordered_map<std::string, double>& mass_by_atom_type();

bool is_standard_protein_residue(const std::string& name);
void configure_xponge_residue_links(ResidueType& residue_type);
void put_template(ResidueType residue_type);
ResidueType residue_type_from_molecule_residue(const Molecule& molecule, const Residue& residue);
void add_minimal_protein_template(const std::string& name);

void upsert_lj_atom_type(const std::string& atom_type, const std::string& lj_type);
void upsert_mass(const std::string& atom_type, double mass);
void upsert_lj_parameter(const std::string& lj_type, double epsilon, double rmin);
void upsert_amber_cmap_key(const std::string& key, std::uint32_t resolution,
                           const std::vector<double>& parameters);

}  // namespace xpongecpp
