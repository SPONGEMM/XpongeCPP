#include "core.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace xpongecpp {
namespace {

std::vector<std::string> split_ws(const std::string& line) {
    std::istringstream input(line);
    std::vector<std::string> words;
    std::string word;
    while (input >> word) {
        words.push_back(word);
    }
    return words;
}

bool is_section_line(const std::vector<std::string>& words) {
    return words.size() >= 2 && !words[1].empty() && words[1][0] == '!';
}

std::vector<int> read_psf_ints(std::istringstream& input, std::size_t count) {
    std::vector<int> values;
    values.reserve(count);
    std::string line;
    while (values.size() < count && std::getline(input, line)) {
        for (const auto& word : split_ws(line)) {
            if (!word.empty() && word[0] == '!') {
                break;
            }
            values.push_back(std::stoi(word));
            if (values.size() == count) {
                break;
            }
        }
    }
    if (values.size() != count) {
        throw std::runtime_error("unexpected EOF while reading PSF integer section");
    }
    return values;
}

std::string residue_key(const std::string& segid, const std::string& resnr, const std::string& resname) {
    return segid + ":" + resnr + ":" + resname;
}

struct PsfAtomSignature {
    std::string type;
    double charge{0.0};
};

using PsfResidueSignature = std::unordered_map<std::string, PsfAtomSignature>;

bool signature_atom_matches(const PsfResidueSignature& signature, const std::string& atom_name,
                            const std::string& atom_type, double charge) {
    const auto it = signature.find(atom_name);
    return it == signature.end() ||
           (it->second.type == atom_type && std::abs(it->second.charge - charge) <= 1e-6);
}

bool signature_matches_residue_atoms(const PsfResidueSignature& signature, const Molecule& molecule,
                                     const Residue& residue) {
    for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
        const auto& atom = molecule.atoms[residue.atom_begin + local];
        if (!signature_atom_matches(signature, atom.name, atom.type, atom.charge)) {
            return false;
        }
    }
    return true;
}

std::string next_residue_type_name(const std::string& base_name,
                                   const std::vector<std::pair<std::string, PsfResidueSignature>>& signatures) {
    for (int suffix = 1;; ++suffix) {
        const std::string candidate = base_name + "_" + std::to_string(suffix);
        const auto found = std::find_if(signatures.begin(), signatures.end(), [&](const auto& entry) {
            return entry.first == candidate;
        });
        if (found == signatures.end()) {
            return candidate;
        }
    }
}

void update_psf_residue_type(Molecule& molecule, Residue& residue, const std::string& atom_name,
                             const std::string& atom_type, double charge,
                             std::unordered_map<std::string, std::vector<std::pair<std::string, PsfResidueSignature>>>&
                                 signatures_by_resname) {
    auto& signatures = signatures_by_resname[residue.name];
    if (signatures.empty()) {
        signatures.push_back({residue.name, {}});
    }

    auto current = std::find_if(signatures.begin(), signatures.end(), [&](const auto& entry) {
        return entry.first == residue.type_name;
    });
    if (current == signatures.end()) {
        current = signatures.begin();
        residue.type_name = current->first;
    }
    if (signature_matches_residue_atoms(current->second, molecule, residue) &&
        signature_atom_matches(current->second, atom_name, atom_type, charge)) {
        current->second[atom_name] = {atom_type, charge};
        return;
    }

    for (auto it = signatures.begin(); it != signatures.end(); ++it) {
        if (signature_matches_residue_atoms(it->second, molecule, residue) &&
            signature_atom_matches(it->second, atom_name, atom_type, charge)) {
            residue.type_name = it->first;
            it->second[atom_name] = {atom_type, charge};
            return;
        }
    }

    PsfResidueSignature new_signature;
    for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
        const auto& atom = molecule.atoms[residue.atom_begin + local];
        new_signature[atom.name] = {atom.type, atom.charge};
    }
    new_signature[atom_name] = {atom_type, charge};
    residue.type_name = next_residue_type_name(residue.name, signatures);
    signatures.push_back({residue.type_name, std::move(new_signature)});
}

void append_psf_atom(Molecule& molecule, std::unordered_map<std::string, ResidueId>& residue_by_key,
                     const std::string& segid, const std::string& resnr, const std::string& resname,
                     const std::string& atom_name, const std::string& atom_type, double charge, double mass,
                     std::unordered_map<std::string, std::vector<std::pair<std::string, PsfResidueSignature>>>&
                         signatures_by_resname) {
    const auto key = residue_key(segid, resnr, resname);
    auto it = residue_by_key.find(key);
    if (it == residue_by_key.end()) {
        const ResidueId residue_id = static_cast<ResidueId>(molecule.residues.size());
        it = residue_by_key.emplace(key, residue_id).first;
        Residue residue;
        residue.name = resname;
        residue.type_name = resname;
        residue.original_name = resname;
        (void)segid;
        residue.pdb_resseq = std::stoi(resnr);
        residue.atom_begin = static_cast<AtomId>(molecule.atoms.size());
        molecule.residues.push_back(std::move(residue));
    }
    auto& residue = molecule.residues[it->second];
    update_psf_residue_type(molecule, residue, atom_name, atom_type, charge, signatures_by_resname);
    Atom atom;
    atom.name = atom_name;
    atom.type = atom_type;
    atom.element = guess_element(atom_name, atom_type);
    atom.residue = it->second;
    atom.charge = charge;
    atom.mass = mass;
    molecule.atoms.push_back(std::move(atom));
    ++residue.atom_count;
}

std::vector<std::vector<AtomId>> adjacency(const Molecule& molecule) {
    std::vector<std::vector<AtomId>> graph(molecule.atoms.size());
    for (const auto& bond : molecule.explicit_bonds) {
        graph[bond.atom1].push_back(bond.atom2);
        graph[bond.atom2].push_back(bond.atom1);
    }
    for (const auto& link : molecule.residue_links) {
        graph[link.atom1].push_back(link.atom2);
        graph[link.atom2].push_back(link.atom1);
    }
    return graph;
}

void split_residues_by_connectivity(Molecule& molecule) {
    const auto graph = adjacency(molecule);
    std::vector<std::vector<AtomId>> residue_components;
    std::vector<bool> linked_atoms(molecule.atoms.size(), false);
    for (const auto& link : molecule.residue_links) {
        linked_atoms[link.atom1] = true;
        linked_atoms[link.atom2] = true;
    }
    for (const auto& residue : molecule.residues) {
        const auto push_whole_residue = [&]() {
            std::vector<AtomId> component;
            component.reserve(residue.atom_count);
            for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
                component.push_back(residue.atom_begin + local);
            }
            residue_components.push_back(std::move(component));
        };
        bool has_linked_atom = false;
        for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
            if (linked_atoms[residue.atom_begin + local]) {
                has_linked_atom = true;
                break;
            }
        }
        if (has_linked_atom || residue.atom_count == 0) {
            push_whole_residue();
            continue;
        }
        std::vector<bool> seen(residue.atom_count, false);
        std::vector<std::vector<AtomId>> components;
        for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
            if (seen[local]) {
                continue;
            }
            std::vector<AtomId> component;
            std::queue<AtomId> queue;
            queue.push(residue.atom_begin + local);
            seen[local] = true;
            while (!queue.empty()) {
                const AtomId atom = queue.front();
                queue.pop();
                component.push_back(atom);
                for (const AtomId next : graph[atom]) {
                    if (molecule.atoms[next].residue != molecule.atoms[atom].residue) {
                        continue;
                    }
                    const std::uint32_t next_local = next - residue.atom_begin;
                    if (!seen[next_local]) {
                        seen[next_local] = true;
                        queue.push(next);
                    }
                }
            }
            components.push_back(std::move(component));
        }
        if (components.size() <= 1) {
            push_whole_residue();
            continue;
        }
        residue_components.insert(residue_components.end(), components.begin(), components.end());
    }

    std::vector<Atom> new_atoms;
    std::vector<Residue> new_residues;
    std::vector<AtomId> old_to_new(molecule.atoms.size(), std::numeric_limits<AtomId>::max());
    new_atoms.reserve(molecule.atoms.size());
    new_residues.reserve(residue_components.size());
    for (const auto& component : residue_components) {
        if (component.empty()) {
            continue;
        }
        Residue residue = molecule.residues[molecule.atoms[component.front()].residue];
        residue.atom_begin = static_cast<AtomId>(new_atoms.size());
        residue.atom_count = static_cast<std::uint32_t>(component.size());
        const ResidueId new_residue_id = static_cast<ResidueId>(new_residues.size());
        for (const AtomId old_atom_id : component) {
            Atom atom = molecule.atoms[old_atom_id];
            atom.residue = new_residue_id;
            old_to_new[old_atom_id] = static_cast<AtomId>(new_atoms.size());
            new_atoms.push_back(std::move(atom));
        }
        new_residues.push_back(std::move(residue));
    }

    for (auto& bond : molecule.explicit_bonds) {
        bond.atom1 = old_to_new.at(bond.atom1);
        bond.atom2 = old_to_new.at(bond.atom2);
    }
    for (auto& link : molecule.residue_links) {
        link.atom1 = old_to_new.at(link.atom1);
        link.atom2 = old_to_new.at(link.atom2);
    }
    molecule.atoms = std::move(new_atoms);
    molecule.residues = std::move(new_residues);
}

std::unordered_map<std::string, Molecule> split_molecules_by_connectivity(const Molecule& molecule) {
    const auto graph = adjacency(molecule);
    std::vector<int> component(molecule.atoms.size(), 0);
    int component_count = 0;
    for (AtomId atom = 0; atom < molecule.atoms.size(); ++atom) {
        if (component[atom] != 0) {
            continue;
        }
        ++component_count;
        std::queue<AtomId> queue;
        queue.push(atom);
        component[atom] = component_count;
        while (!queue.empty()) {
            const AtomId current = queue.front();
            queue.pop();
            for (const AtomId next : graph[current]) {
                if (component[next] == 0) {
                    component[next] = component_count;
                    queue.push(next);
                }
            }
        }
    }

    std::unordered_map<std::string, Molecule> out;
    std::vector<std::unordered_map<AtomId, AtomId>> atom_maps(static_cast<std::size_t>(component_count + 1));
    for (int comp = 1; comp <= component_count; ++comp) {
        out.emplace(molecule.name + "_" + std::to_string(comp), Molecule(molecule.name + "_" + std::to_string(comp)));
    }
    for (const auto& residue : molecule.residues) {
        if (residue.atom_count == 0) {
            continue;
        }
        const int comp = component[residue.atom_begin];
        auto& target = out.at(molecule.name + "_" + std::to_string(comp));
        Residue new_residue = residue;
        new_residue.atom_begin = static_cast<AtomId>(target.atoms.size());
        new_residue.atom_count = 0;
        const ResidueId new_residue_id = static_cast<ResidueId>(target.residues.size());
        for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
            const AtomId old_atom = residue.atom_begin + local;
            if (component[old_atom] != comp) {
                continue;
            }
            Atom atom = molecule.atoms[old_atom];
            atom.residue = new_residue_id;
            atom_maps[comp][old_atom] = static_cast<AtomId>(target.atoms.size());
            target.atoms.push_back(std::move(atom));
            ++new_residue.atom_count;
        }
        target.residues.push_back(std::move(new_residue));
    }
    for (const auto& bond : molecule.explicit_bonds) {
        const int comp = component[bond.atom1];
        if (comp == component[bond.atom2]) {
            auto& map = atom_maps[comp];
            out.at(molecule.name + "_" + std::to_string(comp)).explicit_bonds.push_back({map.at(bond.atom1), map.at(bond.atom2)});
        }
    }
    return out;
}

}  // namespace

PsfData load_molpsf_text(const std::string& text, const std::string& split_by) {
    std::istringstream input(text);
    std::string line;
    if (!std::getline(input, line) || line.find("PSF") == std::string::npos) {
        throw std::runtime_error("Not a PSF file");
    }
    Molecule molecule("psf");
    std::unordered_map<int, AtomId> atom_by_psf_index;
    std::unordered_map<std::string, ResidueId> residue_by_key;
    std::unordered_map<std::string, std::vector<std::pair<std::string, PsfResidueSignature>>> signatures_by_resname;

    while (std::getline(input, line)) {
        const auto words = split_ws(line);
        if (!is_section_line(words)) {
            continue;
        }
        const int count = std::stoi(words[0]);
        const std::string section = words[1];
        if (section.find("NTITLE") != std::string::npos) {
            for (int i = 0; i < count; ++i) {
                std::getline(input, line);
            }
        } else if (section.find("NATOM") != std::string::npos) {
            for (int i = 0; i < count; ++i) {
                if (!std::getline(input, line)) {
                    throw std::runtime_error("unexpected EOF in PSF NATOM section");
                }
                const auto atom_words = split_ws(line);
                if (atom_words.size() < 8) {
                    throw std::runtime_error("bad PSF atom line");
                }
                append_psf_atom(molecule, residue_by_key, atom_words[1], atom_words[2], atom_words[3],
                                atom_words[4], atom_words[5], std::stod(atom_words[6]), std::stod(atom_words[7]),
                                signatures_by_resname);
                atom_by_psf_index[std::stoi(atom_words[0])] = static_cast<AtomId>(molecule.atoms.size() - 1);
            }
        } else if (section.find("NBOND") != std::string::npos) {
            const auto ints = read_psf_ints(input, static_cast<std::size_t>(count * 2));
            for (std::size_t i = 0; i + 1 < ints.size(); i += 2) {
                const AtomId atom1 = atom_by_psf_index.at(ints[i]);
                const AtomId atom2 = atom_by_psf_index.at(ints[i + 1]);
                molecule.explicit_bonds.push_back({atom1, atom2});
                if (molecule.atoms[atom1].residue != molecule.atoms[atom2].residue) {
                    molecule.add_residue_link(atom1, atom2);
                }
            }
        }
    }
    split_residues_by_connectivity(molecule);
    PsfData data{molecule, {}};
    if (split_by == "connectivity") {
        data.molecules = split_molecules_by_connectivity(data.molecule);
    } else {
        data.molecules.emplace(data.molecule.name, data.molecule);
    }
    return data;
}

}  // namespace xpongecpp
