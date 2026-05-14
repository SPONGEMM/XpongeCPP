#pragma once

#include "core.hpp"

#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace xpongecpp {

struct AtomTypeInfo {
    double mass{0.0};
    double epsilon{0.0};
    double rmin{0.0};
};

struct CharmmUreyParameter {
    std::array<std::string, 3> types;
    double k{0.0};
    double theta{0.0};
    double k_ub{0.0};
    double r13{0.0};
};

struct NBfixParameter {
    std::string type1;
    std::string type2;
    double a{0.0};
    double b{0.0};
    double kee{0.0};
};

std::unordered_map<std::string, AtomTypeInfo>& atom_type_registry();
std::vector<CharmmUreyParameter>& charmm_urey_registry();
std::vector<NBfixParameter>& nbfix_registry();

std::string trim_copy(const std::string& input);
std::string strip_comment(const std::string& line);
std::vector<std::string> split_ws(const std::string& line);
std::string upper_copy(std::string text);
AtomId one_based_atom(const std::string& token);

void register_external_lj(const std::string& type, double epsilon, double rmin);
void register_external_mass(const std::string& type, double mass);

void append_atom(Molecule& molecule, std::unordered_map<int, ResidueId>& residue_id_by_number,
                 int residue_number, const std::string& residue_name, const std::string& atom_name,
                 const std::string& atom_type, double charge, double mass);

}  // namespace xpongecpp
