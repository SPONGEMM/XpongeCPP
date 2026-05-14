#include "pdb_internal.hpp"
#include "pdb_records.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace xpongecpp {

std::string pdb_trimmed_copy(const std::string& input) {
    return pdb_trim_copy(input);
}

std::string pdb_upper_copy(std::string value) {
    for (auto& ch : value) {
        ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    }
    return value;
}

char char_at(const std::string& line, std::size_t pos, char fallback) {
    return pos < line.size() ? line[pos] : fallback;
}

std::string pdb_string(const std::string& line, std::size_t pos, std::size_t len) {
    if (line.size() <= pos) {
        return "";
    }
    return pdb_trimmed_copy(line.substr(pos, std::min(len, line.size() - pos)));
}

double pdb_float(const std::string& line, std::size_t pos, std::size_t len, const char* field_name,
                 double default_value) {
    if (line.size() <= pos) {
        if (!std::isnan(default_value)) {
            return default_value;
        }
        throw std::invalid_argument(std::string("missing PDB coordinate field: ") + field_name);
    }
    const auto field = pdb_trimmed_copy(line.substr(pos, std::min(len, line.size() - pos)));
    if (field.empty()) {
        if (!std::isnan(default_value)) {
            return default_value;
        }
        throw std::invalid_argument(std::string("empty PDB coordinate field: ") + field_name);
    }
    try {
        return std::stod(field);
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string("invalid PDB coordinate field ") + field_name + ": " + field);
    }
}

std::unordered_map<std::string, std::string>& pdb_head_map() {
    static std::unordered_map<std::string, std::string> map;
    return map;
}

std::unordered_map<std::string, std::string>& pdb_tail_map() {
    static std::unordered_map<std::string, std::string> map;
    return map;
}

std::unordered_map<std::string, std::string>& pdb_save_map() {
    static std::unordered_map<std::string, std::string> map;
    return map;
}

std::unordered_map<std::string, std::string>& pdb_alias_map() {
    static std::unordered_map<std::string, std::string> map;
    return map;
}

std::unordered_map<std::string, HisNames>& his_map() {
    static std::unordered_map<std::string, HisNames> map;
    return map;
}

std::unordered_set<std::string> protein_residue_names() {
    return {"ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN", "GLU", "GLY", "HID", "HIE",
            "HIP", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
}

bool is_hydrogen_name(const std::string& atom_name) {
    return !atom_name.empty() &&
           (atom_name[0] == 'H' ||
            (atom_name.size() > 1 && (atom_name[0] == '1' || atom_name[0] == '2' || atom_name[0] == '3') &&
             atom_name[1] == 'H'));
}

void set_atom_defaults(Atom& atom, const std::string& residue_name) {
    const auto residue_upper = pdb_upper_copy(residue_name);
    const auto atom_upper = pdb_upper_copy(atom.name);
    if (residue_upper == "WAT" || residue_upper == "HOH" || residue_upper == "TIP3") {
        atom.type = atom.element == "H" ? "HW" : "OW";
        atom.charge = atom.element == "H" ? 0.417 : -0.834;
    } else if (residue_upper == "NA" || residue_upper == "Na+") {
        atom.type = "Na+";
        atom.element = "Na";
        atom.charge = 1.0;
    } else if (residue_upper == "CL" || residue_upper == "Cl-") {
        atom.type = "Cl-";
        atom.element = "Cl";
        atom.charge = -1.0;
    } else if (atom_upper == "OXT") {
        atom.type = "O";
    } else {
        atom.type = atom.element;
    }
    atom.mass = default_mass_for_element(atom.element);
}

int hy36_decode(int width, const std::string& field) {
    return pdb_hy36_decode(width, field);
}

std::optional<int> pdb_int(const std::string& line, std::size_t pos, std::size_t len) {
    if (line.size() <= pos) {
        return std::nullopt;
    }
    const auto field = line.substr(pos, std::min(len, line.size() - pos));
    const auto text = pdb_trimmed_copy(field);
    if (text.empty()) {
        return std::nullopt;
    }
    try {
        return hy36_decode(static_cast<int>(len), field);
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

std::string pdb_int_field(int width, int value) {
    return pdb_hy36_field(width, value);
}

std::vector<std::string> pdb_split_ws(const std::string& line) {
    std::istringstream input(line);
    std::vector<std::string> out;
    std::string word;
    while (input >> word) {
        out.push_back(word);
    }
    return out;
}

std::optional<std::tuple<char, int, char>> ssbond_ref(const std::string& line, bool second) {
    const char chain = second ? char_at(line, 29) : char_at(line, 15);
    const auto resseq = second ? pdb_int(line, 31, 4) : pdb_int(line, 17, 4);
    if (resseq && chain != ' ') {
        return std::tuple<char, int, char>{chain, *resseq, ' '};
    }
    const auto words = pdb_split_ws(line);
    const std::size_t chain_index = second ? 6 : 3;
    const std::size_t resseq_index = second ? 7 : 4;
    if (words.size() > resseq_index) {
        return std::tuple<char, int, char>{words[chain_index].empty() ? ' ' : words[chain_index][0],
                                           std::stoi(words[resseq_index]), ' '};
    }
    return std::nullopt;
}

ResidueSelectorSets parse_unterminal_residues(const std::vector<std::string>& selectors) {
    ResidueSelectorSets out;
    for (auto selector : selectors) {
        selector = pdb_trimmed_copy(selector);
        if (selector.empty()) {
            continue;
        }
        if (selector.front() == '(' || selector.front() == '[') {
            selector.erase(std::remove_if(selector.begin(), selector.end(),
                                          [](char c) { return c == '(' || c == ')' || c == '[' || c == ']' ||
                                                             c == '\'' || c == '"' || c == ' '; }),
                           selector.end());
            std::replace(selector.begin(), selector.end(), ',', ':');
        }
        char chain = '\0';
        std::string residue_text = selector;
        const auto colon = selector.find(':');
        if (colon != std::string::npos) {
            chain = selector[0];
            residue_text = selector.substr(colon + 1);
        }
        char insertion = ' ';
        if (residue_text.size() > 1 && std::isalpha(static_cast<unsigned char>(residue_text.back()))) {
            insertion = residue_text.back();
            residue_text.pop_back();
        }
        const int resseq = std::stoi(residue_text);
        if (chain == '\0') {
            out.all_resseq.insert(resseq);
        } else if (insertion == ' ') {
            out.chain_resseq.insert({chain, resseq});
        } else {
            out.chain_resseq_ins.insert({chain, resseq, insertion});
        }
    }
    return out;
}

bool is_unterminal(const ResidueSelectorSets& selectors, char chain_id, int resseq, char insertion_code) {
    return selectors.all_resseq.count(resseq) != 0 ||
           selectors.chain_resseq.count({chain_id, resseq}) != 0 ||
           selectors.chain_resseq_ins.count({chain_id, resseq, insertion_code}) != 0;
}

void apply_template_atom_properties(Molecule& molecule, bool ignore_unknown_name) {
    std::vector<Atom> new_atoms;
    new_atoms.reserve(molecule.atoms.size());
    for (auto& residue : molecule.residues) {
        const AtomId new_begin = static_cast<AtomId>(new_atoms.size());
        std::uint32_t new_count = 0;
        const ResidueType* residue_type = nullptr;
        if (has_template(residue.name)) {
            residue_type = &get_residue_template(residue.name);
        }
        for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
            auto atom = molecule.atoms[residue.atom_begin + local];
            if (residue_type) {
                try {
                    const auto template_index = residue_type->atom_index(atom.name);
                    const auto& template_atom = residue_type->atoms()[template_index];
                    atom.type = template_atom.type;
                    atom.element = template_atom.element;
                    atom.charge = template_atom.charge;
                    atom.mass = template_atom.mass;
                } catch (const std::exception&) {
                    if (ignore_unknown_name) {
                        continue;
                    }
                    throw;
                }
            }
            atom.residue = static_cast<ResidueId>(&residue - molecule.residues.data());
            new_atoms.push_back(std::move(atom));
            ++new_count;
        }
        residue.atom_begin = new_begin;
        residue.atom_count = new_count;
    }
    molecule.atoms = std::move(new_atoms);
}

AtomId find_atom(const Molecule& molecule, const Residue& residue, const std::string& name) {
    for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
        const AtomId atom_id = residue.atom_begin + local;
        if (molecule.atoms[atom_id].name == name) {
            return atom_id;
        }
    }
    return static_cast<AtomId>(molecule.atoms.size());
}

void add_template_residue_link(Molecule& molecule, ResidueId left_id, ResidueId right_id) {
    const auto& left = molecule.residues[left_id];
    const auto& right = molecule.residues[right_id];
    if (!has_template(left.name) || !has_template(right.name)) {
        return;
    }
    const auto& left_type = get_residue_template(left.name);
    const auto& right_type = get_residue_template(right.name);
    if (left_type.tail().empty() || right_type.head().empty()) {
        return;
    }
    const AtomId tail = find_atom(molecule, left, left_type.tail());
    const AtomId head = find_atom(molecule, right, right_type.head());
    if (tail < molecule.atoms.size() && head < molecule.atoms.size()) {
        molecule.add_residue_link(tail, head);
    }
}

std::string residue_key(std::uint32_t segment, char chain, int resseq, char insertion, const std::string& name) {
    return std::to_string(segment) + ":" + chain + ":" + std::to_string(resseq) + ":" + insertion + ":" + name;
}

std::string normalized_residue_name(const std::string& name) {
    const auto it = pdb_save_map().find(name);
    std::string out = it == pdb_save_map().end() ? name : it->second;
    if (out.size() > 3 && out.find("__") != std::string::npos) {
        out = out.substr(0, out.find("__"));
    }
    return out.size() > 3 ? out.substr(0, 3) : out;
}

bool atom_pair_equals(const Molecule& molecule, const ResidueLink& link, const Residue& a, const std::string& atom_a,
                      const Residue& b, const std::string& atom_b) {
    const AtomId id_a = find_atom(molecule, a, atom_a);
    const AtomId id_b = find_atom(molecule, b, atom_b);
    return id_a < molecule.atoms.size() && id_b < molecule.atoms.size() &&
           ((link.atom1 == id_a && link.atom2 == id_b) || (link.atom1 == id_b && link.atom2 == id_a));
}

void register_pdb_residue_name_mapping(const std::string& place, const std::string& pdb_name,
                                       const std::string& real_name) {
    if (place == "head") {
        pdb_head_map()[pdb_name] = real_name;
    } else if (place == "tail") {
        pdb_tail_map()[pdb_name] = real_name;
    } else {
        throw std::invalid_argument("PDB residue name mapping place should be head or tail");
    }
    pdb_save_map()[real_name] = pdb_name;
}

void register_pdb_residue_alias_mapping(const std::string& pdb_name, const std::string& real_name) {
    pdb_alias_map()[pdb_name] = real_name;
}

void register_his_mapping(const std::string& residue_name, const std::string& hid, const std::string& hie,
                          const std::string& hip) {
    his_map()[residue_name] = {hid, hie, hip};
}

}  // namespace xpongecpp
