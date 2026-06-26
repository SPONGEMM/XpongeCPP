#include "pdb_internal.hpp"

#include <algorithm>
#include <limits>
#include <map>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace xpongecpp {

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

}  // namespace xpongecpp
