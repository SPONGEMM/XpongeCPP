#include "core.hpp"
#include "pdb_records.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

namespace xpongecpp {
namespace {

std::string trim_copy(const std::string& input) {
    return pdb_trim_copy(input);
}

std::string upper_copy(std::string value) {
    for (auto& ch : value) {
        ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    }
    return value;
}

char char_at(const std::string& line, std::size_t pos, char fallback = ' ') {
    return pos < line.size() ? line[pos] : fallback;
}

std::string pdb_string(const std::string& line, std::size_t pos, std::size_t len) {
    if (line.size() <= pos) {
        return "";
    }
    return trim_copy(line.substr(pos, std::min(len, line.size() - pos)));
}

double pdb_float(const std::string& line, std::size_t pos, std::size_t len, const char* field_name,
                 double default_value = std::numeric_limits<double>::quiet_NaN()) {
    if (line.size() <= pos) {
        if (!std::isnan(default_value)) {
            return default_value;
        }
        throw std::invalid_argument(std::string("missing PDB coordinate field: ") + field_name);
    }
    const auto field = trim_copy(line.substr(pos, std::min(len, line.size() - pos)));
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

struct HisNames {
    std::string hid;
    std::string hie;
    std::string hip;
};

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
    const auto residue_upper = upper_copy(residue_name);
    const auto atom_upper = upper_copy(atom.name);
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
    const auto text = trim_copy(field);
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

struct ResidueSelectorSets {
    std::unordered_set<int> all_resseq;
    std::set<std::pair<char, int>> chain_resseq;
    std::set<std::tuple<char, int, char>> chain_resseq_ins;
};

std::vector<std::string> split_ws(const std::string& line) {
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
    const auto words = split_ws(line);
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
        selector = trim_copy(selector);
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

}  // namespace

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

Molecule load_pdb_text(const std::string& text) {
    return load_pdb_text(text, PdbLoadOptions{});
}

Molecule load_pdb_text(const std::string& text, const PdbLoadOptions& options) {
    Molecule molecule("PDB");
    std::istringstream input(text);
    std::string line;
    std::string current_residue_key;
    ResidueId current_residue_id = 0;
    std::uint32_t current_segment = 0;
    bool segment_has_residue = false;
    std::unordered_set<char> processed_chains;
    char last_chain = ' ';
    int insertion_count = 0;
    std::vector<std::string> ssbond_lines;
    std::vector<std::string> link_lines;
    std::vector<std::string> conect_lines;
    std::unordered_map<int, AtomId> serial_to_atom;
    std::map<std::tuple<char, int, char>, ResidueId> chain_residue_to_id;
    std::vector<bool> residue_unterminal;
    std::unordered_set<ResidueId> oxt_residues;
    const auto unterminal = parse_unterminal_residues(options.unterminal_residues);

    const auto finish_segment = [&]() {
        if (!molecule.residues.empty()) {
            auto& last = molecule.residues.back();
            if (!residue_unterminal.back()) {
                const auto it = pdb_tail_map().find(last.name);
                if (it != pdb_tail_map().end()) {
                    last.name = it->second;
                    last.type_name = last.name;
                }
            }
        }
        processed_chains.insert(last_chain);
        ++current_segment;
        segment_has_residue = false;
        current_residue_key.clear();
        insertion_count = 0;
    };

    while (std::getline(input, line)) {
        if (options.read_cryst1 && line.rfind("CRYST1", 0) == 0) {
            try {
                molecule.box_length = {std::stod(trim_copy(line.substr(6, 9))),
                                       std::stod(trim_copy(line.substr(15, 9))),
                                       std::stod(trim_copy(line.substr(24, 9)))};
                molecule.box_angle = {std::stod(trim_copy(line.substr(33, 7))),
                                      std::stod(trim_copy(line.substr(40, 7))),
                                      std::stod(trim_copy(line.substr(47, 7)))};
                molecule.has_box = true;
            } catch (const std::exception&) {
                throw std::invalid_argument("invalid CRYST1 record");
            }
            continue;
        }
        if (line.rfind("SSBOND", 0) == 0) {
            ssbond_lines.push_back(line);
            continue;
        }
        if (line.rfind("LINK", 0) == 0) {
            link_lines.push_back(line);
            continue;
        }
        if (!options.ignore_conect && line.rfind("CONECT", 0) == 0) {
            conect_lines.push_back(line);
            continue;
        }
        if (line.rfind("TER", 0) == 0) {
            finish_segment();
            continue;
        }
        if (!(line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0)) {
            continue;
        }

        auto resseq_opt = pdb_int(line, 22, 4);
        auto serial_opt = pdb_int(line, 6, 5);
        if (!resseq_opt || !serial_opt) {
            continue;
        }
        int resseq = *resseq_opt + insertion_count;
        const char insertion = char_at(line, 26);
        char chain = char_at(line, 21);
        const char chain_in_file = chain == ' ' ? ' ' : chain;
        if (chain == ' ' || processed_chains.count(chain) != 0) {
            chain = ' ';
        }
        last_chain = char_at(line, 21);
        std::string residue_name = upper_copy(pdb_string(line, 17, 3));
        if (const auto alias = pdb_alias_map().find(residue_name); alias != pdb_alias_map().end()) {
            residue_name = alias->second;
        }
        const std::string atom_name = pdb_string(line, 12, 4);
        const char altloc = char_at(line, 16);
        if (altloc != ' ' && altloc != options.position_need) {
            continue;
        }
        if (options.ignore_hydrogen && is_hydrogen_name(atom_name)) {
            continue;
        }
        const bool skip_terminal = is_unterminal(unterminal, chain_in_file, *resseq_opt, insertion);
        const auto key = residue_key(current_segment, chain, resseq, insertion, residue_name);
        if (molecule.residues.empty() || key != current_residue_key) {
            if (!current_residue_key.empty() && insertion != ' ') {
                ++insertion_count;
                ++resseq;
            }
            current_residue_key = residue_key(current_segment, chain, resseq, insertion, residue_name);
            current_residue_id = static_cast<ResidueId>(molecule.residues.size());
            Residue residue;
            residue.name = residue_name.empty() ? "UNK" : residue_name;
            residue.original_name = residue.name;
            if (!segment_has_residue && !skip_terminal && pdb_head_map().count(residue.name) != 0) {
                residue.name = pdb_head_map().at(residue.name);
            }
            residue.type_name = residue.name;
            residue.chain_id = chain_in_file;
            residue.effective_chain_id = chain;
            residue.segment_id = current_segment;
            residue.pdb_resseq = *resseq_opt;
            residue.insertion_code = insertion;
            residue.is_hetero = line.rfind("HETATM", 0) == 0;
            residue.atom_begin = static_cast<AtomId>(molecule.atoms.size());
            molecule.residues.push_back(residue);
            segment_has_residue = true;
            residue_unterminal.push_back(skip_terminal);
            chain_residue_to_id[{chain_in_file, *resseq_opt, insertion}] = current_residue_id;
        }

        Atom atom;
        atom.name = atom_name;
        atom.serial = *serial_opt;
        atom.altloc = altloc;
        atom.occupancy = pdb_float(line, 54, 6, "occupancy", 1.0);
        atom.temp_factor = pdb_float(line, 60, 6, "temperature factor", 0.0);
        atom.record_name = line.rfind("HETATM", 0) == 0 ? "HETATM" : "ATOM";
        atom.element = guess_element(atom_name, pdb_string(line, 76, 2));
        atom.residue = current_residue_id;
        atom.x = pdb_float(line, 30, 8, "x");
        atom.y = pdb_float(line, 38, 8, "y");
        atom.z = pdb_float(line, 46, 8, "z");
        set_atom_defaults(atom, molecule.residues[current_residue_id].name);
        serial_to_atom[atom.serial] = static_cast<AtomId>(molecule.atoms.size());
        molecule.atoms.push_back(std::move(atom));
        molecule.residues[current_residue_id].atom_count += 1;
        if (atom_name == "OXT") {
            oxt_residues.insert(current_residue_id);
        }
    }
    if (!molecule.residues.empty() && !current_residue_key.empty()) {
        finish_segment();
    }

    for (ResidueId residue_id = 0; residue_id < molecule.residues.size(); ++residue_id) {
        auto& residue = molecule.residues[residue_id];
        if (options.judge_histone) {
            const auto it = his_map().find(residue.name);
            if (it != his_map().end()) {
                bool delta = false;
                bool epsilon = false;
                for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
                    const auto& atom = molecule.atoms[residue.atom_begin + local];
                    delta = delta || atom.name == "HD1";
                    epsilon = epsilon || atom.name == "HE2";
                }
                residue.name = delta ? (epsilon ? it->second.hip : it->second.hid) : it->second.hie;
                residue.type_name = residue.name;
            }
        }
        if (oxt_residues.count(residue_id) != 0 && !residue_unterminal[residue_id]) {
            const auto proteins = protein_residue_names();
            if (proteins.count(residue.name) != 0 && pdb_tail_map().count(residue.name) != 0) {
                residue.name = pdb_tail_map().at(residue.name);
                residue.type_name = residue.name;
            } else if (residue.name.size() > 1 && residue.name[0] == 'N' &&
                       proteins.count(residue.name.substr(1)) != 0 && pdb_tail_map().count(residue.name.substr(1)) != 0) {
                residue.name = "C" + residue.name.substr(1);
                residue.type_name = residue.name;
            }
        }
    }

    const auto apply_ssbond_name = [&](ResidueId residue_id) {
        auto& residue = molecule.residues[residue_id];
        const auto& head = pdb_head_map();
        const auto& tail = pdb_tail_map();
        const bool is_head = std::any_of(head.begin(), head.end(),
                                         [&](const auto& item) { return item.second == residue.name; });
        const bool is_tail = std::any_of(tail.begin(), tail.end(),
                                         [&](const auto& item) { return item.second == residue.name; });
        residue.name = is_head ? "NCYX" : (is_tail ? "CCYX" : "CYX");
        residue.type_name = residue.name;
    };

    for (const auto& ssbond : ssbond_lines) {
        const auto ref_a = ssbond_ref(ssbond, false);
        const auto ref_b = ssbond_ref(ssbond, true);
        if (!ref_a || !ref_b) {
            continue;
        }
        const auto it_a = chain_residue_to_id.find(*ref_a);
        const auto it_b = chain_residue_to_id.find(*ref_b);
        if (it_a == chain_residue_to_id.end() || it_b == chain_residue_to_id.end()) {
            continue;
        }
        apply_ssbond_name(it_a->second);
        apply_ssbond_name(it_b->second);
    }

    apply_template_atom_properties(molecule, options.ignore_unknown_name);

    for (ResidueId residue_id = 0; residue_id + 1 < molecule.residues.size(); ++residue_id) {
        if (molecule.residues[residue_id].segment_id == molecule.residues[residue_id + 1].segment_id) {
            add_template_residue_link(molecule, residue_id, residue_id + 1);
        }
    }

    for (const auto& ssbond : ssbond_lines) {
        const auto ref_a = ssbond_ref(ssbond, false);
        const auto ref_b = ssbond_ref(ssbond, true);
        if (!ref_a || !ref_b) {
            continue;
        }
        const auto it_a = chain_residue_to_id.find(*ref_a);
        const auto it_b = chain_residue_to_id.find(*ref_b);
        if (it_a == chain_residue_to_id.end() || it_b == chain_residue_to_id.end()) {
            continue;
        }
        const auto& residue_a = molecule.residues[it_a->second];
        const auto& residue_b = molecule.residues[it_b->second];
        const auto& type_a = get_residue_template(residue_a.name);
        const auto& type_b = get_residue_template(residue_b.name);
        const auto atom_name_a = type_a.connect_atoms().at("ssbond");
        const auto atom_name_b = type_b.connect_atoms().at("ssbond");
        const AtomId atom_a = find_atom(molecule, residue_a, atom_name_a);
        const AtomId atom_b = find_atom(molecule, residue_b, atom_name_b);
        if (atom_a < molecule.atoms.size() && atom_b < molecule.atoms.size()) {
            molecule.add_residue_link(atom_a, atom_b);
        }
    }

    for (const auto& link : link_lines) {
        const std::string atom_name_a = pdb_string(link, 12, 4);
        const std::string atom_name_b = pdb_string(link, 42, 4);
        const char chain_a = char_at(link, 21);
        const char chain_b = char_at(link, 51);
        const auto res_a = pdb_int(link, 22, 4);
        const auto res_b = pdb_int(link, 52, 4);
        if (!res_a || !res_b) {
            continue;
        }
        const auto it_a = chain_residue_to_id.find({chain_a, *res_a, char_at(link, 26)});
        const auto it_b = chain_residue_to_id.find({chain_b, *res_b, char_at(link, 56)});
        if (it_a == chain_residue_to_id.end() || it_b == chain_residue_to_id.end()) {
            continue;
        }
        const AtomId atom_a = find_atom(molecule, molecule.residues[it_a->second], atom_name_a);
        const AtomId atom_b = find_atom(molecule, molecule.residues[it_b->second], atom_name_b);
        if (atom_a < molecule.atoms.size() && atom_b < molecule.atoms.size()) {
            molecule.add_residue_link(atom_a, atom_b);
        }
    }

    for (const auto& connect : conect_lines) {
        const auto atom = pdb_int(connect, 6, 5);
        if (!atom || serial_to_atom.count(*atom) == 0) {
            continue;
        }
        for (std::size_t pos = 11; pos + 5 <= connect.size(); pos += 5) {
            const auto other = pdb_int(connect, pos, 5);
            if (!other || serial_to_atom.count(*other) == 0) {
                continue;
            }
            const AtomId atom_a = serial_to_atom.at(*atom);
            const AtomId atom_b = serial_to_atom.at(*other);
            if (molecule.atoms[atom_a].residue != molecule.atoms[atom_b].residue) {
                molecule.add_residue_link(atom_a, atom_b);
            }
        }
    }

    return molecule;
}

void save_pdb(const Molecule& molecule, const std::filesystem::path& filename) {
    if (!molecule.validate()) {
        throw std::invalid_argument("cannot export invalid molecule");
    }
    std::vector<char> chain_ids(molecule.residues.size(), ' ');
    std::map<char, std::map<int, int>> chains;
    char alphabet = 'A';
    std::size_t start = 0;
    const auto residue_link_exists = [&](ResidueId left_id, ResidueId right_id) {
        const auto& left = molecule.residues[left_id];
        const auto& right = molecule.residues[right_id];
        if (!has_template(left.name) || !has_template(right.name)) {
            return false;
        }
        const auto& left_type = get_residue_template(left.name);
        const auto& right_type = get_residue_template(right.name);
        if (left_type.tail().empty() || right_type.head().empty()) {
            return false;
        }
        for (const auto& link : molecule.residue_links) {
            if (atom_pair_equals(molecule, link, left, left_type.tail(), right, right_type.head())) {
                return true;
            }
        }
        return false;
    };
    for (std::size_t i = 1; i <= molecule.residues.size(); ++i) {
        bool new_chain = i == molecule.residues.size();
        if (!new_chain) {
            new_chain = !residue_link_exists(static_cast<ResidueId>(i - 1), static_cast<ResidueId>(i));
        }
        if (new_chain) {
            const auto length = i - start;
            if (length > 1) {
                for (std::size_t j = 0; j < length; ++j) {
                    chain_ids[start + j] = alphabet;
                    chains[alphabet][static_cast<int>(j + 1)] = static_cast<int>(start + j + 1);
                }
                ++alphabet;
            }
            start = i;
        }
    }

    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("failed to open PDB output: " + filename.string());
    }
    out << "REMARK   Generated By Xponge (Molecule)\n";
    out << std::fixed << std::setprecision(3);
    if (molecule.has_box) {
        out << "CRYST1" << std::setw(9) << molecule.box_length[0] << std::setw(9) << molecule.box_length[1]
            << std::setw(9) << molecule.box_length[2] << std::setw(7) << molecule.box_angle[0] << std::setw(7)
            << molecule.box_angle[1] << std::setw(7) << molecule.box_angle[2] << " P 1           1\n";
    }
    for (const auto& [chain_id, index_map] : chains) {
        std::vector<std::string> names;
        for (const auto& [pdb_index, residue_index_1] : index_map) {
            (void)pdb_index;
            names.push_back(normalized_residue_name(molecule.residues[residue_index_1 - 1].name));
        }
        for (std::size_t row = 0; row < names.size(); row += 13) {
            out << "SEQRES " << std::setw(3) << (row / 13 + 1) << " " << chain_id << std::setw(5)
                << names.size() << "  ";
            for (std::size_t i = row; i < std::min(row + 13, names.size()); ++i) {
                out << std::setw(3) << names[i] << (i + 1 == std::min(row + 13, names.size()) ? "" : " ");
            }
            out << "\n";
        }
    }

    std::vector<int> residue_pdb_indices(molecule.residues.size(), 1);
    std::unordered_map<char, int> next_index_by_chain;
    for (std::size_t i = 0; i < molecule.residues.size(); ++i) {
        if (chain_ids[i] == ' ') {
            residue_pdb_indices[i] = 1;
        } else {
            residue_pdb_indices[i] = ++next_index_by_chain[chain_ids[i]];
        }
    }

    struct LinkRecord {
        char chain_a{' '};
        int resseq_a{0};
        char chain_b{' '};
        int resseq_b{0};
        std::string atom_a;
        std::string resname_a;
        std::string atom_b;
        std::string resname_b;
    };
    std::map<AtomId, std::vector<AtomId>> connects;
    std::vector<std::tuple<char, int, char, int>> ssbonds;
    std::vector<LinkRecord> links;
    const auto residue_atom_name = [&](const Residue& residue, AtomId atom_id) {
        if (atom_id < residue.atom_begin || atom_id >= residue.atom_begin + residue.atom_count) {
            return std::string{};
        }
        return molecule.atoms[atom_id].name;
    };
    const auto is_ssbond_atom = [&](const Residue& residue, AtomId atom_id) {
        if (!has_template(residue.name)) {
            return false;
        }
        const auto& type = get_residue_template(residue.name);
        const auto it = type.connect_atoms().find("ssbond");
        return it != type.connect_atoms().end() && residue_atom_name(residue, atom_id) == it->second;
    };
    for (const auto& link : molecule.residue_links) {
        AtomId atom_a = link.atom1;
        AtomId atom_b = link.atom2;
        if (atom_a > atom_b) {
            std::swap(atom_a, atom_b);
        }
        const auto res_a = molecule.atoms[atom_a].residue;
        const auto res_b = molecule.atoms[atom_b].residue;
        const char chain_a = chain_ids[res_a];
        const char chain_b = chain_ids[res_b];
        if (res_a + 1 == res_b && chain_a == chain_b && chain_a != ' ') {
            continue;
        }
        if (chain_a == ' ' || chain_b == ' ') {
            connects[atom_a].push_back(atom_b);
            connects[atom_b].push_back(atom_a);
            continue;
        }
        const auto& residue_a = molecule.residues[res_a];
        const auto& residue_b = molecule.residues[res_b];
        const std::string resname_a = normalized_residue_name(residue_a.name);
        const std::string resname_b = normalized_residue_name(residue_b.name);
        if (resname_a == "CYX" && resname_b == "CYX" &&
            is_ssbond_atom(residue_a, atom_a) && is_ssbond_atom(residue_b, atom_b)) {
            ssbonds.push_back({chain_a, residue_pdb_indices[res_a], chain_b, residue_pdb_indices[res_b]});
        } else {
            links.push_back({chain_a, residue_pdb_indices[res_a], chain_b, residue_pdb_indices[res_b],
                             residue_atom_name(residue_a, atom_a), resname_a,
                             residue_atom_name(residue_b, atom_b), resname_b});
        }
    }
    std::sort(ssbonds.begin(), ssbonds.end());
    for (std::size_t i = 0; i < ssbonds.size(); ++i) {
        const auto& [chain_a, resseq_a, chain_b, resseq_b] = ssbonds[i];
        out << "SSBOND " << std::setw(3) << (i + 1) << " CYX " << chain_a << pdb_int_field(4, resseq_a)
            << "    CYX " << chain_b << pdb_int_field(4, resseq_b) << "\n";
    }
    std::sort(links.begin(), links.end(), [](const LinkRecord& lhs, const LinkRecord& rhs) {
        return std::tie(lhs.chain_a, lhs.resseq_a, lhs.chain_b, lhs.resseq_b, lhs.atom_a, lhs.atom_b) <
               std::tie(rhs.chain_a, rhs.resseq_a, rhs.chain_b, rhs.resseq_b, rhs.atom_a, rhs.atom_b);
    });
    for (const auto& link : links) {
        out << "LINK        " << std::setw(4) << std::left << link.atom_a.substr(0, 4) << std::right << " "
            << std::setw(3) << link.resname_a << " " << link.chain_a << pdb_int_field(4, link.resseq_a)
            << "                " << std::setw(4) << std::left << link.atom_b.substr(0, 4) << std::right << " "
            << std::setw(3) << link.resname_b << " " << link.chain_b << pdb_int_field(4, link.resseq_b) << "\n";
    }

    std::size_t serial = 1;
    for (std::size_t residue_index = 0; residue_index < molecule.residues.size(); ++residue_index) {
        const auto& residue = molecule.residues[residue_index];
        const int resseq = residue_pdb_indices[residue_index];
        for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
            const auto& atom = molecule.atoms[residue.atom_begin + local];
            const auto resname = normalized_residue_name(residue.name);
            out << (atom.record_name == "HETATM" ? "HETATM" : "ATOM  ") << pdb_int_field(5, static_cast<int>(serial))
                << " " << std::left << std::setw(4) << atom.name.substr(0, 4) << std::right << " " << std::setw(3)
                << resname << " " << chain_ids[residue_index] << pdb_int_field(4, resseq) << "    "
                << std::setw(8) << atom.x << std::setw(8) << atom.y << std::setw(8) << atom.z << std::setw(6)
                << std::setprecision(2) << atom.occupancy << std::setw(6) << atom.temp_factor << std::setprecision(3)
                << "          " << std::setw(2) << atom.element.substr(0, 2) << "\n";
            ++serial;
        }
        const bool is_last = residue_index + 1 == molecule.residues.size();
        if (is_last || chain_ids[residue_index + 1] != chain_ids[residue_index] || chain_ids[residue_index] == ' ') {
            out << "TER\n";
        }
    }

    for (auto& [atom, atoms] : connects) {
        std::sort(atoms.begin(), atoms.end());
        for (std::size_t i = 0; i < atoms.size(); i += 4) {
            out << "CONECT" << pdb_int_field(5, static_cast<int>(atom + 1));
            for (std::size_t j = i; j < std::min(i + 4, atoms.size()); ++j) {
                out << pdb_int_field(5, static_cast<int>(atoms[j] + 1));
            }
            out << "\n";
        }
    }
    out << "END\n";
}

}  // namespace xpongecpp
