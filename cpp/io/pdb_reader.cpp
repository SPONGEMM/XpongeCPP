#include "pdb_internal.hpp"

#include <algorithm>
#include <cctype>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace xpongecpp {
namespace {

std::string mmcif_lower_copy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

std::string mmcif_clean(const std::string& value) {
    const auto out = pdb_trimmed_copy(value);
    return out == "." || out == "?" ? std::string{} : out;
}

std::optional<int> mmcif_int(const std::string& value) {
    const auto cleaned = mmcif_clean(value);
    if (cleaned.empty()) {
        return std::nullopt;
    }
    try {
        return std::stoi(cleaned);
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

std::optional<double> mmcif_double(const std::string& value) {
    const auto cleaned = mmcif_clean(value);
    if (cleaned.empty()) {
        return std::nullopt;
    }
    try {
        return std::stod(cleaned);
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

char mmcif_char(const std::string& value, char fallback = ' ') {
    const auto cleaned = mmcif_clean(value);
    return cleaned.empty() ? fallback : cleaned[0];
}

std::vector<std::string> mmcif_tokenize(const std::string& text) {
    std::vector<std::string> tokens;
    std::istringstream input(text);
    std::string line;
    while (std::getline(input, line)) {
        if (!line.empty() && line.front() == ';') {
            std::string block = line.substr(1);
            bool closed = false;
            while (std::getline(input, line)) {
                if (!line.empty() && line.front() == ';') {
                    closed = true;
                    break;
                }
                block += "\n";
                block += line;
            }
            if (!closed) {
                throw std::invalid_argument("unterminated mmCIF text field");
            }
            tokens.push_back(std::move(block));
            continue;
        }
        std::string current;
        char quote = '\0';
        for (std::size_t i = 0; i < line.size(); ++i) {
            const char c = line[i];
            if (quote != '\0') {
                if (c == quote && (i + 1 == line.size() || std::isspace(static_cast<unsigned char>(line[i + 1])))) {
                    tokens.push_back(current);
                    current.clear();
                    quote = '\0';
                } else {
                    current.push_back(c);
                }
                continue;
            }
            if (c == '#') {
                break;
            }
            if (std::isspace(static_cast<unsigned char>(c))) {
                if (!current.empty()) {
                    tokens.push_back(current);
                    current.clear();
                }
                continue;
            }
            if (c == '\'' || c == '"') {
                if (current.empty()) {
                    quote = c;
                } else {
                    current.push_back(c);
                }
                continue;
            }
            current.push_back(c);
        }
        if (quote != '\0') {
            throw std::invalid_argument("unterminated mmCIF quoted value");
        }
        if (!current.empty()) {
            tokens.push_back(std::move(current));
        }
    }
    return tokens;
}

using MmcifData = std::unordered_map<std::string, std::vector<std::string>>;
using MmcifRow = std::unordered_map<std::string, std::string>;

MmcifData mmcif_parse(const std::string& text) {
    const auto tokens = mmcif_tokenize(text);
    MmcifData data;
    std::size_t i = 0;
    while (i < tokens.size()) {
        const auto lower = mmcif_lower_copy(tokens[i]);
        if (lower.rfind("data_", 0) == 0 || lower.rfind("save_", 0) == 0) {
            ++i;
            continue;
        }
        if (lower == "loop_") {
            ++i;
            std::vector<std::string> tags;
            while (i < tokens.size() && !tokens[i].empty() && tokens[i][0] == '_') {
                tags.push_back(mmcif_lower_copy(tokens[i]));
                ++i;
            }
            if (tags.empty()) {
                throw std::invalid_argument("mmCIF loop without tags");
            }
            std::vector<std::string> values;
            while (i < tokens.size()) {
                const auto next_lower = mmcif_lower_copy(tokens[i]);
                const bool row_boundary = values.size() % tags.size() == 0;
                if (row_boundary &&
                    (next_lower == "loop_" || next_lower.rfind("data_", 0) == 0 ||
                     next_lower.rfind("save_", 0) == 0 || (!tokens[i].empty() && tokens[i][0] == '_'))) {
                    break;
                }
                values.push_back(tokens[i]);
                ++i;
            }
            if (values.size() % tags.size() != 0) {
                throw std::invalid_argument("mmCIF loop row has incomplete values");
            }
            for (std::size_t tag_index = 0; tag_index < tags.size(); ++tag_index) {
                auto& column = data[tags[tag_index]];
                for (std::size_t value_index = tag_index; value_index < values.size(); value_index += tags.size()) {
                    column.push_back(values[value_index]);
                }
            }
            continue;
        }
        if (!tokens[i].empty() && tokens[i][0] == '_') {
            if (i + 1 >= tokens.size()) {
                throw std::invalid_argument("mmCIF tag without value: " + tokens[i]);
            }
            data[mmcif_lower_copy(tokens[i])] = {tokens[i + 1]};
            i += 2;
            continue;
        }
        ++i;
    }
    return data;
}

std::vector<MmcifRow> mmcif_rows(const MmcifData& data, const std::string& category) {
    const std::string prefix = "_" + mmcif_lower_copy(category) + ".";
    std::vector<std::string> tags;
    std::size_t count = 0;
    for (const auto& [tag, values] : data) {
        if (tag.rfind(prefix, 0) == 0) {
            tags.push_back(tag);
            count = std::max(count, values.size());
        }
    }
    std::vector<MmcifRow> rows;
    rows.reserve(count);
    for (std::size_t row_index = 0; row_index < count; ++row_index) {
        MmcifRow row;
        for (const auto& tag : tags) {
            const auto& values = data.at(tag);
            row[tag] = row_index < values.size() ? values[row_index] : std::string{};
        }
        rows.push_back(std::move(row));
    }
    return rows;
}

std::string mmcif_first(const MmcifRow& row, std::initializer_list<const char*> tags) {
    for (const auto* tag : tags) {
        const auto it = row.find(mmcif_lower_copy(tag));
        if (it == row.end()) {
            continue;
        }
        const auto value = mmcif_clean(it->second);
        if (!value.empty()) {
            return value;
        }
    }
    return {};
}

std::string mmcif_value(const MmcifRow& row, const std::string& tag) {
    const auto it = row.find(mmcif_lower_copy(tag));
    return it == row.end() ? std::string{} : mmcif_clean(it->second);
}

std::string mmcif_scalar(const MmcifData& data, const std::string& tag) {
    const auto it = data.find(mmcif_lower_copy(tag));
    if (it == data.end() || it->second.empty()) {
        return {};
    }
    return mmcif_clean(it->second.front());
}

std::string mmcif_atom_key(const std::string& asym, const std::string& seq, const std::string& comp,
                           const std::string& atom, char insertion) {
    return mmcif_clean(asym) + "|" + mmcif_clean(seq) + "|" + pdb_upper_copy(mmcif_clean(comp)) + "|" +
           pdb_upper_copy(mmcif_clean(atom)) + "|" + std::string(1, insertion);
}

struct MmcifAtomInfo {
    std::string site_id;
    std::string atom_name;
    std::string element;
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double occupancy{1.0};
    double temp_factor{0.0};
    std::string record_name{"ATOM"};
    char altloc{' '};
    std::string label_key;
    std::string auth_key;
};

struct MmcifResidueInfo {
    char chain_id{' '};
    int resseq{0};
    char insertion_code{' '};
    std::string base_name;
    std::string name;
    std::string label_asym;
    std::string label_seq;
    std::string label_comp;
    std::vector<MmcifAtomInfo> atoms;
    std::unordered_set<std::string> atom_names;
    bool has_oxt{false};
    bool disulfide{false};
};

void mmcif_add_explicit_bond(Molecule& molecule, AtomId atom1, AtomId atom2) {
    if (atom1 >= molecule.atoms.size() || atom2 >= molecule.atoms.size() || atom1 == atom2) {
        return;
    }
    const auto lo = std::min(atom1, atom2);
    const auto hi = std::max(atom1, atom2);
    for (const auto& bond : molecule.explicit_bonds) {
        if (std::min(bond.atom1, bond.atom2) == lo && std::max(bond.atom1, bond.atom2) == hi) {
            return;
        }
    }
    molecule.explicit_bonds.push_back({lo, hi});
}

void mmcif_add_connection(Molecule& molecule, AtomId atom1, AtomId atom2) {
    if (atom1 >= molecule.atoms.size() || atom2 >= molecule.atoms.size() || atom1 == atom2) {
        return;
    }
    if (molecule.atoms[atom1].residue == molecule.atoms[atom2].residue) {
        mmcif_add_explicit_bond(molecule, atom1, atom2);
        return;
    }
    molecule.add_residue_link(atom1, atom2);
}

void mmcif_remove_connection(Molecule& molecule, AtomId atom1, AtomId atom2) {
    if (atom1 >= molecule.atoms.size() || atom2 >= molecule.atoms.size() || atom1 == atom2) {
        return;
    }
    const auto lo = std::min(atom1, atom2);
    const auto hi = std::max(atom1, atom2);
    if (molecule.atoms[atom1].residue == molecule.atoms[atom2].residue) {
        molecule.explicit_bonds.erase(
            std::remove_if(molecule.explicit_bonds.begin(), molecule.explicit_bonds.end(),
                           [&](const auto& bond) {
                               return std::min(bond.atom1, bond.atom2) == lo &&
                                      std::max(bond.atom1, bond.atom2) == hi;
                           }),
            molecule.explicit_bonds.end());
        return;
    }
    molecule.residue_links.erase(
        std::remove_if(molecule.residue_links.begin(), molecule.residue_links.end(),
                       [&](const auto& link) {
                           return std::min(link.atom1, link.atom2) == lo &&
                                  std::max(link.atom1, link.atom2) == hi;
                       }),
        molecule.residue_links.end());
}

std::string mmcif_histidine_name(const std::string& residue_name, const std::unordered_set<std::string>& atom_names) {
    const auto it = his_map().find(residue_name);
    if (it == his_map().end()) {
        return residue_name;
    }
    const bool delta = atom_names.count("HD1") != 0;
    const bool epsilon = atom_names.count("HE2") != 0;
    return delta ? (epsilon ? it->second.hip : it->second.hid) : it->second.hie;
}

}  // namespace

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
    std::vector<std::pair<bool, bool>> residue_terminal;
    std::unordered_set<ResidueId> oxt_residues;
    const auto unterminal = parse_unterminal_residues(options.unterminal_residues);

    const auto finish_segment = [&]() {
        if (!molecule.residues.empty()) {
            auto& last = molecule.residues.back();
            if (options.infer_terminals && !residue_unterminal.back()) {
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
                molecule.box_length = {std::stod(pdb_trimmed_copy(line.substr(6, 9))),
                                       std::stod(pdb_trimmed_copy(line.substr(15, 9))),
                                       std::stod(pdb_trimmed_copy(line.substr(24, 9)))};
                molecule.box_angle = {std::stod(pdb_trimmed_copy(line.substr(33, 7))),
                                      std::stod(pdb_trimmed_copy(line.substr(40, 7))),
                                      std::stod(pdb_trimmed_copy(line.substr(47, 7)))};
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
        std::string residue_name = pdb_upper_copy(pdb_string(line, 17, 3));
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
        const auto explicit_terminal =
            terminal_residue_flags(options.terminal_residues, chain_in_file, *resseq_opt, insertion);
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
            const bool should_map_head =
                explicit_terminal.first || (options.infer_terminals && !segment_has_residue && !skip_terminal);
            if (should_map_head && pdb_head_map().count(residue.name) != 0) {
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
            residue_terminal.push_back(explicit_terminal);
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
        atom.x = pdb_float(line, 30, 8, "x", std::numeric_limits<double>::quiet_NaN());
        atom.y = pdb_float(line, 38, 8, "y", std::numeric_limits<double>::quiet_NaN());
        atom.z = pdb_float(line, 46, 8, "z", std::numeric_limits<double>::quiet_NaN());
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
        if (options.infer_terminals && oxt_residues.count(residue_id) != 0 && !residue_unterminal[residue_id]) {
            residue.name = tail_mapped_residue_name(residue.name);
            residue.type_name = residue.name;
        }
        if (residue_terminal[residue_id].second) {
            residue.name = tail_mapped_residue_name(residue.name);
            residue.type_name = residue.name;
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

Molecule load_mmcif_text(const std::string& text) {
    return load_mmcif_text(text, MmcifLoadOptions{});
}

Molecule load_mmcif_text(const std::string& text, const MmcifLoadOptions& options) {
    const auto data = mmcif_parse(text);
    const auto atom_rows = mmcif_rows(data, "atom_site");
    if (atom_rows.empty()) {
        throw std::invalid_argument("mmCIF file does not contain _atom_site records");
    }

    std::unordered_set<std::string> model_values;
    for (const auto& row : atom_rows) {
        auto model = mmcif_first(row, {"_atom_site.pdbx_pdb_model_num"});
        model_values.insert(model.empty() ? "1" : model);
    }
    if (!options.model_id && model_values.size() > 1) {
        throw std::invalid_argument("mmCIF contains multiple models; pass model_id explicitly");
    }
    const auto selected_model = options.model_id;
    const auto unterminal = parse_unterminal_residues(options.unterminal_residues);

    std::vector<MmcifResidueInfo> residues;
    std::string current_key;
    for (const auto& row : atom_rows) {
        auto row_model = mmcif_first(row, {"_atom_site.pdbx_pdb_model_num"});
        if (row_model.empty()) {
            row_model = "1";
        }
        if (selected_model && row_model != *selected_model) {
            continue;
        }
        const auto altloc_text = mmcif_first(row, {"_atom_site.label_alt_id"});
        const char altloc = altloc_text.empty() ? ' ' : altloc_text[0];
        if (altloc != ' ' && altloc != options.position_need) {
            continue;
        }
        const std::string atom_name =
            mmcif_first(row, {"_atom_site.auth_atom_id", "_atom_site.label_atom_id"});
        if (options.ignore_hydrogen && is_hydrogen_name(atom_name)) {
            continue;
        }
        std::string residue_name =
            pdb_upper_copy(mmcif_first(row, {"_atom_site.auth_comp_id", "_atom_site.label_comp_id"}));
        if (const auto alias = pdb_alias_map().find(residue_name); alias != pdb_alias_map().end()) {
            residue_name = alias->second;
        }
        const std::string auth_seq = mmcif_first(row, {"_atom_site.auth_seq_id"});
        const std::string label_seq = mmcif_first(row, {"_atom_site.label_seq_id"});
        auto resseq = auth_seq.empty() ? std::optional<int>{} : mmcif_int(auth_seq);
        if (!resseq) {
            resseq = mmcif_int(label_seq);
        }
        if (!resseq) {
            resseq = static_cast<int>(residues.size()) + 1;
        }
        const std::string auth_asym = mmcif_first(row, {"_atom_site.auth_asym_id", "_atom_site.label_asym_id"});
        const char chain_id = mmcif_char(auth_asym);
        const char insertion = mmcif_char(mmcif_first(row, {"_atom_site.pdbx_pdb_ins_code"}));
        const auto key = residue_key(0, chain_id, *resseq, insertion, residue_name);
        if (key != current_key) {
            current_key = key;
            MmcifResidueInfo info;
            info.chain_id = chain_id;
            info.resseq = *resseq;
            info.insertion_code = insertion;
            info.base_name = residue_name;
            info.name = residue_name;
            info.label_asym = mmcif_first(row, {"_atom_site.label_asym_id"});
            info.label_seq = label_seq;
            info.label_comp = mmcif_first(row, {"_atom_site.label_comp_id"});
            residues.push_back(std::move(info));
        }
        MmcifAtomInfo atom_info;
        atom_info.site_id = mmcif_first(row, {"_atom_site.id"});
        atom_info.atom_name = atom_name;
        atom_info.element = mmcif_first(row, {"_atom_site.type_symbol"});
        const auto x = mmcif_double(mmcif_first(row, {"_atom_site.cartn_x"}));
        const auto y = mmcif_double(mmcif_first(row, {"_atom_site.cartn_y"}));
        const auto z = mmcif_double(mmcif_first(row, {"_atom_site.cartn_z"}));
        if (!x || !y || !z) {
            throw std::invalid_argument("mmCIF atom_site row has invalid Cartesian coordinates");
        }
        atom_info.x = *x;
        atom_info.y = *y;
        atom_info.z = *z;
        atom_info.occupancy = mmcif_double(mmcif_first(row, {"_atom_site.occupancy"})).value_or(1.0);
        atom_info.temp_factor = mmcif_double(mmcif_first(row, {"_atom_site.b_iso_or_equiv"})).value_or(0.0);
        atom_info.record_name = mmcif_first(row, {"_atom_site.group_pdb"});
        if (atom_info.record_name.empty()) {
            atom_info.record_name = "ATOM";
        }
        atom_info.altloc = altloc;
        atom_info.label_key = mmcif_atom_key(
            mmcif_first(row, {"_atom_site.label_asym_id"}),
            label_seq,
            mmcif_first(row, {"_atom_site.label_comp_id"}),
            mmcif_first(row, {"_atom_site.label_atom_id"}),
            insertion);
        atom_info.auth_key = mmcif_atom_key(
            auth_asym,
            auth_seq.empty() ? label_seq : auth_seq,
            mmcif_first(row, {"_atom_site.auth_comp_id", "_atom_site.label_comp_id"}),
            atom_name,
            insertion);
        residues.back().atom_names.insert(atom_name);
        if (atom_name == "OXT") {
            residues.back().has_oxt = true;
        }
        residues.back().atoms.push_back(std::move(atom_info));
    }
    if (residues.empty()) {
        throw std::invalid_argument("mmCIF model selection produced no atoms");
    }

    std::unordered_map<std::string, ResidueId> label_residue_to_id;
    std::unordered_map<std::string, ResidueId> auth_residue_to_id;
    for (ResidueId i = 0; i < residues.size(); ++i) {
        label_residue_to_id[mmcif_atom_key(residues[i].label_asym, residues[i].label_seq,
                                           residues[i].label_comp, "", residues[i].insertion_code)] = i;
        auth_residue_to_id[mmcif_atom_key(std::string(1, residues[i].chain_id), std::to_string(residues[i].resseq),
                                          residues[i].base_name, "", residues[i].insertion_code)] = i;
    }
    for (const auto& row : mmcif_rows(data, "struct_conn")) {
        const auto atom1 = pdb_upper_copy(
            mmcif_first(row, {"_struct_conn.ptnr1_auth_atom_id", "_struct_conn.ptnr1_label_atom_id"}));
        const auto atom2 = pdb_upper_copy(
            mmcif_first(row, {"_struct_conn.ptnr2_auth_atom_id", "_struct_conn.ptnr2_label_atom_id"}));
        const auto comp1 = pdb_upper_copy(
            mmcif_first(row, {"_struct_conn.ptnr1_auth_comp_id", "_struct_conn.ptnr1_label_comp_id"}));
        const auto comp2 = pdb_upper_copy(
            mmcif_first(row, {"_struct_conn.ptnr2_auth_comp_id", "_struct_conn.ptnr2_label_comp_id"}));
        if (atom1 != "SG" || atom2 != "SG" || comp1 != "CYS" || comp2 != "CYS") {
            continue;
        }
        const auto label_key1 = mmcif_atom_key(
            mmcif_first(row, {"_struct_conn.ptnr1_label_asym_id"}),
            mmcif_first(row, {"_struct_conn.ptnr1_label_seq_id"}),
            mmcif_first(row, {"_struct_conn.ptnr1_label_comp_id"}),
            "",
            mmcif_char(mmcif_first(row, {"_struct_conn.pdbx_ptnr1_pdb_ins_code"})));
        const auto label_key2 = mmcif_atom_key(
            mmcif_first(row, {"_struct_conn.ptnr2_label_asym_id"}),
            mmcif_first(row, {"_struct_conn.ptnr2_label_seq_id"}),
            mmcif_first(row, {"_struct_conn.ptnr2_label_comp_id"}),
            "",
            mmcif_char(mmcif_first(row, {"_struct_conn.pdbx_ptnr2_pdb_ins_code"})));
        const auto auth_key1 = mmcif_atom_key(
            mmcif_first(row, {"_struct_conn.ptnr1_auth_asym_id", "_struct_conn.ptnr1_label_asym_id"}),
            mmcif_first(row, {"_struct_conn.ptnr1_auth_seq_id", "_struct_conn.ptnr1_label_seq_id"}),
            mmcif_first(row, {"_struct_conn.ptnr1_auth_comp_id", "_struct_conn.ptnr1_label_comp_id"}),
            "",
            mmcif_char(mmcif_first(row, {"_struct_conn.pdbx_ptnr1_pdb_ins_code"})));
        const auto auth_key2 = mmcif_atom_key(
            mmcif_first(row, {"_struct_conn.ptnr2_auth_asym_id", "_struct_conn.ptnr2_label_asym_id"}),
            mmcif_first(row, {"_struct_conn.ptnr2_auth_seq_id", "_struct_conn.ptnr2_label_seq_id"}),
            mmcif_first(row, {"_struct_conn.ptnr2_auth_comp_id", "_struct_conn.ptnr2_label_comp_id"}),
            "",
            mmcif_char(mmcif_first(row, {"_struct_conn.pdbx_ptnr2_pdb_ins_code"})));
        const auto mark_disulfide = [&](const std::string& auth_key, const std::string& label_key) {
            if (const auto it = auth_residue_to_id.find(auth_key); it != auth_residue_to_id.end()) {
                residues[it->second].disulfide = true;
            } else if (const auto fallback = label_residue_to_id.find(label_key); fallback != label_residue_to_id.end()) {
                residues[fallback->second].disulfide = true;
            }
        };
        mark_disulfide(auth_key1, label_key1);
        mark_disulfide(auth_key2, label_key2);
    }

    std::unordered_map<char, ResidueId> chain_first;
    std::unordered_map<char, ResidueId> chain_last;
    for (ResidueId i = 0; i < residues.size(); ++i) {
        chain_first.emplace(residues[i].chain_id, i);
        chain_last[residues[i].chain_id] = i;
    }
    for (ResidueId i = 0; i < residues.size(); ++i) {
        auto& info = residues[i];
        std::string name = options.judge_histone ? mmcif_histidine_name(info.base_name, info.atom_names)
                                                 : info.base_name;
        const bool skip_terminal = is_unterminal(unterminal, info.chain_id, info.resseq, info.insertion_code);
        const auto explicit_terminal =
            terminal_residue_flags(options.terminal_residues, info.chain_id, info.resseq, info.insertion_code);
        const bool should_map_head =
            explicit_terminal.first ||
            (options.infer_terminals && !skip_terminal && chain_first.at(info.chain_id) == i);
        const bool should_map_tail =
            explicit_terminal.second ||
            (options.infer_terminals && !skip_terminal &&
             (chain_last.at(info.chain_id) == i || info.has_oxt));
        if (should_map_head && pdb_head_map().count(name) != 0) {
            name = pdb_head_map().at(name);
        }
        if (should_map_tail) {
            name = tail_mapped_residue_name(name);
        }
        if (info.disulfide) {
            const auto& head = pdb_head_map();
            const auto& tail = pdb_tail_map();
            const bool is_head = std::any_of(head.begin(), head.end(),
                                             [&](const auto& item) { return item.second == name; });
            const bool is_tail = std::any_of(tail.begin(), tail.end(),
                                             [&](const auto& item) { return item.second == name; });
            name = is_head ? "NCYX" : (is_tail ? "CCYX" : "CYX");
        }
        info.name = name;
    }

    Molecule molecule("mmCIF");
    for (ResidueId residue_id = 0; residue_id < residues.size(); ++residue_id) {
        const auto& info = residues[residue_id];
        Residue residue;
        residue.name = info.name.empty() ? "UNK" : info.name;
        residue.type_name = residue.name;
        residue.original_name = info.base_name;
        residue.chain_id = info.chain_id;
        residue.effective_chain_id = info.chain_id;
        residue.segment_id = static_cast<std::uint32_t>(info.chain_id);
        residue.pdb_resseq = info.resseq;
        residue.insertion_code = info.insertion_code;
        residue.is_hetero = false;
        residue.atom_begin = static_cast<AtomId>(molecule.atoms.size());
        for (const auto& atom_info : info.atoms) {
            Atom atom;
            atom.name = atom_info.atom_name;
            atom.serial = mmcif_int(atom_info.site_id).value_or(static_cast<int>(molecule.atoms.size()) + 1);
            atom.altloc = atom_info.altloc;
            atom.occupancy = atom_info.occupancy;
            atom.temp_factor = atom_info.temp_factor;
            atom.record_name = atom_info.record_name == "HETATM" ? "HETATM" : "ATOM";
            atom.element = guess_element(atom.name, atom_info.element);
            atom.residue = residue_id;
            atom.x = atom_info.x;
            atom.y = atom_info.y;
            atom.z = atom_info.z;
            set_atom_defaults(atom, residue.name);
            molecule.atoms.push_back(std::move(atom));
            ++residue.atom_count;
        }
        molecule.residues.push_back(std::move(residue));
    }

    apply_template_atom_properties(molecule, options.ignore_unknown_name);

    std::unordered_map<std::string, AtomId> atom_by_site;
    std::unordered_map<std::string, AtomId> atom_by_label;
    std::unordered_map<std::string, AtomId> atom_by_auth;
    std::unordered_map<std::string, AtomId> atom_by_external;
    for (ResidueId residue_id = 0; residue_id < residues.size(); ++residue_id) {
        const auto& residue = molecule.residues[residue_id];
        const auto& info = residues[residue_id];
        for (const auto& atom_info : info.atoms) {
            const AtomId atom_id = find_atom(molecule, residue, atom_info.atom_name);
            if (atom_id >= molecule.atoms.size()) {
                continue;
            }
            if (!atom_info.site_id.empty()) {
                atom_by_site[atom_info.site_id] = atom_id;
            }
            atom_by_label[atom_info.label_key] = atom_id;
            atom_by_auth[atom_info.auth_key] = atom_id;
            atom_by_external[std::string(1, info.chain_id) + "|" + std::to_string(info.resseq) + "|" +
                             std::string(1, info.insertion_code) + "|" + pdb_upper_copy(info.base_name) + "|" +
                             pdb_upper_copy(molecule.atoms[atom_id].name)] = atom_id;
        }
    }

    for (ResidueId residue_id = 1; residue_id < molecule.residues.size(); ++residue_id) {
        if (molecule.residues[residue_id - 1].chain_id == molecule.residues[residue_id].chain_id) {
            add_template_residue_link(molecule, residue_id - 1, residue_id);
        }
    }

    std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> chem_comp_bonds;
    for (const auto& row : mmcif_rows(data, "chem_comp_bond")) {
        const auto comp = pdb_upper_copy(mmcif_first(row, {"_chem_comp_bond.comp_id"}));
        const auto atom1 = pdb_upper_copy(mmcif_first(row, {"_chem_comp_bond.atom_id_1"}));
        const auto atom2 = pdb_upper_copy(mmcif_first(row, {"_chem_comp_bond.atom_id_2"}));
        if (!comp.empty() && !atom1.empty() && !atom2.empty()) {
            chem_comp_bonds[comp].push_back({atom1, atom2});
        }
    }
    for (ResidueId residue_id = 0; residue_id < molecule.residues.size(); ++residue_id) {
        const auto& info = residues[residue_id];
        const auto& residue = molecule.residues[residue_id];
        const auto it = chem_comp_bonds.find(pdb_upper_copy(info.base_name));
        if (it == chem_comp_bonds.end()) {
            continue;
        }
        for (const auto& [atom1_name, atom2_name] : it->second) {
            const AtomId atom1 = find_atom(molecule, residue, atom1_name);
            const AtomId atom2 = find_atom(molecule, residue, atom2_name);
            mmcif_add_connection(molecule, atom1, atom2);
        }
    }

    const auto resolve_struct_conn_atom = [&](const MmcifRow& row, const std::string& partner) {
        const auto insertion = mmcif_char(mmcif_value(row, "_struct_conn.pdbx_" + partner + "_pdb_ins_code"));
        const auto label_key = mmcif_atom_key(
            mmcif_value(row, "_struct_conn." + partner + "_label_asym_id"),
            mmcif_value(row, "_struct_conn." + partner + "_label_seq_id"),
            mmcif_value(row, "_struct_conn." + partner + "_label_comp_id"),
            mmcif_value(row, "_struct_conn." + partner + "_label_atom_id"),
            insertion);
        auto auth_asym = mmcif_value(row, "_struct_conn." + partner + "_auth_asym_id");
        auto auth_seq = mmcif_value(row, "_struct_conn." + partner + "_auth_seq_id");
        auto auth_comp = mmcif_value(row, "_struct_conn." + partner + "_auth_comp_id");
        auto auth_atom = mmcif_value(row, "_struct_conn." + partner + "_auth_atom_id");
        if (auth_asym.empty()) auth_asym = mmcif_value(row, "_struct_conn." + partner + "_label_asym_id");
        if (auth_seq.empty()) auth_seq = mmcif_value(row, "_struct_conn." + partner + "_label_seq_id");
        if (auth_comp.empty()) auth_comp = mmcif_value(row, "_struct_conn." + partner + "_label_comp_id");
        if (auth_atom.empty()) auth_atom = mmcif_value(row, "_struct_conn." + partner + "_label_atom_id");
        const auto auth_key = mmcif_atom_key(auth_asym, auth_seq, auth_comp, auth_atom, insertion);
        const auto auth_it = atom_by_auth.find(auth_key);
        if (auth_it != atom_by_auth.end()) {
            return auth_it->second;
        }
        const auto label_it = atom_by_label.find(label_key);
        return label_it == atom_by_label.end() ? static_cast<AtomId>(molecule.atoms.size()) : label_it->second;
    };
    for (const auto& row : mmcif_rows(data, "struct_conn")) {
        mmcif_add_connection(molecule, resolve_struct_conn_atom(row, "ptnr1"), resolve_struct_conn_atom(row, "ptnr2"));
    }
    for (const auto& row : mmcif_rows(data, "mokda_bond_semantic")) {
        const auto it1 = atom_by_site.find(mmcif_first(row, {"_mokda_bond_semantic.atom_site_id_1"}));
        const auto it2 = atom_by_site.find(mmcif_first(row, {"_mokda_bond_semantic.atom_site_id_2"}));
        if (it1 != atom_by_site.end() && it2 != atom_by_site.end()) {
            mmcif_add_connection(molecule, it1->second, it2->second);
        }
    }
    for (const auto& row : mmcif_rows(data, "mokda_edit_operation")) {
        const auto op = mmcif_lower_copy(mmcif_first(row, {"_mokda_edit_operation.operation_type"}));
        const auto it1 = atom_by_site.find(mmcif_first(row, {"_mokda_edit_operation.atom_site_id_1"}));
        const auto it2 = atom_by_site.find(mmcif_first(row, {"_mokda_edit_operation.atom_site_id_2"}));
        if (it1 == atom_by_site.end() || it2 == atom_by_site.end()) {
            continue;
        }
        if (op == "delete_bond") {
            mmcif_remove_connection(molecule, it1->second, it2->second);
        } else if (op == "add_bond" || op == "create_bond" || op == "update_bond") {
            mmcif_add_connection(molecule, it1->second, it2->second);
        }
    }

    for (const auto& link : options.residue_links) {
        const auto key1 = std::string(1, link.atom1.chain_id) + "|" + std::to_string(link.atom1.resseq) + "|" +
                          std::string(1, link.atom1.insertion_code) + "|" +
                          pdb_upper_copy(link.atom1.residue_name) + "|" + pdb_upper_copy(link.atom1.atom_name);
        const auto key2 = std::string(1, link.atom2.chain_id) + "|" + std::to_string(link.atom2.resseq) + "|" +
                          std::string(1, link.atom2.insertion_code) + "|" +
                          pdb_upper_copy(link.atom2.residue_name) + "|" + pdb_upper_copy(link.atom2.atom_name);
        const auto it1 = atom_by_external.find(key1);
        const auto it2 = atom_by_external.find(key2);
        if (it1 == atom_by_external.end() || it2 == atom_by_external.end()) {
            throw std::invalid_argument("external residue link selector did not match mmCIF atoms");
        }
        mmcif_add_connection(molecule, it1->second, it2->second);
    }

    if (options.read_cell) {
        const auto a = mmcif_double(mmcif_scalar(data, "_cell.length_a"));
        const auto b = mmcif_double(mmcif_scalar(data, "_cell.length_b"));
        const auto c = mmcif_double(mmcif_scalar(data, "_cell.length_c"));
        const auto alpha = mmcif_double(mmcif_scalar(data, "_cell.angle_alpha"));
        const auto beta = mmcif_double(mmcif_scalar(data, "_cell.angle_beta"));
        const auto gamma = mmcif_double(mmcif_scalar(data, "_cell.angle_gamma"));
        if (a && b && c && alpha && beta && gamma) {
            molecule.box_length = {*a, *b, *c};
            molecule.box_angle = {*alpha, *beta, *gamma};
            molecule.has_box = true;
        }
    }
    return molecule;
}

}  // namespace xpongecpp
