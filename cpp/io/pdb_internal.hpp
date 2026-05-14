#pragma once

#include "core.hpp"

#include <cstddef>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace xpongecpp {

std::string pdb_trimmed_copy(const std::string& input);
std::string pdb_upper_copy(std::string value);
char char_at(const std::string& line, std::size_t pos, char fallback = ' ');
std::string pdb_string(const std::string& line, std::size_t pos, std::size_t len);
double pdb_float(const std::string& line, std::size_t pos, std::size_t len, const char* field_name,
                 double default_value);
int hy36_decode(int width, const std::string& field);
std::optional<int> pdb_int(const std::string& line, std::size_t pos, std::size_t len);
std::string pdb_int_field(int width, int value);

std::unordered_map<std::string, std::string>& pdb_head_map();
std::unordered_map<std::string, std::string>& pdb_tail_map();
std::unordered_map<std::string, std::string>& pdb_save_map();
std::unordered_map<std::string, std::string>& pdb_alias_map();

struct HisNames {
    std::string hid;
    std::string hie;
    std::string hip;
};

std::unordered_map<std::string, HisNames>& his_map();
std::unordered_set<std::string> protein_residue_names();
bool is_hydrogen_name(const std::string& atom_name);
void set_atom_defaults(Atom& atom, const std::string& residue_name);

struct ResidueSelectorSets {
    std::unordered_set<int> all_resseq;
    std::set<std::pair<char, int>> chain_resseq;
    std::set<std::tuple<char, int, char>> chain_resseq_ins;
};

std::vector<std::string> pdb_split_ws(const std::string& line);
std::optional<std::tuple<char, int, char>> ssbond_ref(const std::string& line, bool second);
ResidueSelectorSets parse_unterminal_residues(const std::vector<std::string>& selectors);
bool is_unterminal(const ResidueSelectorSets& selectors, char chain_id, int resseq, char insertion_code);
void apply_template_atom_properties(Molecule& molecule, bool ignore_unknown_name);
AtomId find_atom(const Molecule& molecule, const Residue& residue, const std::string& name);
void add_template_residue_link(Molecule& molecule, ResidueId left_id, ResidueId right_id);
std::string residue_key(std::uint32_t segment, char chain, int resseq, char insertion, const std::string& name);
std::string normalized_residue_name(const std::string& name);
bool atom_pair_equals(const Molecule& molecule, const ResidueLink& link, const Residue& a, const std::string& atom_a,
                      const Residue& b, const std::string& atom_b);

}  // namespace xpongecpp
