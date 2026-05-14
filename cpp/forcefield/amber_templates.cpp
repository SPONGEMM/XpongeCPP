#include "amber_internal.hpp"

#include <algorithm>
#include <fstream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace xpongecpp {

ResidueType residue_type_from_molecule_residue(const Molecule& molecule, const Residue& residue) {
    ResidueType residue_type(residue.name);
    for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
        const auto& atom = molecule.atoms[residue.atom_begin + local];
        residue_type.add_atom(atom.name, atom.type, atom.x, atom.y, atom.z, atom.charge, atom.mass);
    }
    for (const auto& bond : molecule.explicit_bonds) {
        if (bond.atom1 < residue.atom_begin || bond.atom1 >= residue.atom_begin + residue.atom_count ||
            bond.atom2 < residue.atom_begin || bond.atom2 >= residue.atom_begin + residue.atom_count) {
            continue;
        }
        residue_type.add_connectivity(molecule.atoms[bond.atom1].name, molecule.atoms[bond.atom2].name);
    }
    return residue_type;
}

void register_residue_templates_from_mol2_text(const std::string& text) {
    const auto molecule = load_mol2_text(text);
    std::vector<std::pair<ResidueId, ResidueType>> residue_types;
    residue_types.reserve(molecule.residues.size());
    std::unordered_map<ResidueId, std::size_t> residue_to_type;

    for (ResidueId residue_id = 0; residue_id < molecule.residues.size(); ++residue_id) {
        const auto& residue = molecule.residues[residue_id];
        residue_to_type[residue_id] = residue_types.size();
        residue_types.emplace_back(residue_id, ResidueType(residue.name));
        auto& residue_type = residue_types.back().second;
        for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
            const auto& atom = molecule.atoms[residue.atom_begin + local];
            residue_type.add_atom(atom.name, atom.type, atom.x, atom.y, atom.z, atom.charge, atom.mass);
        }
    }

    for (const auto& bond : molecule.explicit_bonds) {
        const auto res1 = molecule.atoms[bond.atom1].residue;
        const auto res2 = molecule.atoms[bond.atom2].residue;
        if (res1 != res2) {
            continue;
        }
        auto& residue_type = residue_types[residue_to_type.at(res1)].second;
        try {
            residue_type.add_connectivity(molecule.atoms[bond.atom1].name, molecule.atoms[bond.atom2].name);
        } catch (const std::exception&) {
        }
    }

    for (auto& residue_type : residue_types) {
        put_template(std::move(residue_type.second));
    }
}

void register_residue_templates_from_mol2_file(const std::filesystem::path& filename) {
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open mol2 template file: " + filename.string());
    }
    std::ostringstream buffer;
    buffer << input.rdbuf();
    register_residue_templates_from_mol2_text(buffer.str());
}

void register_template_molecule_from_mol2_file(const std::filesystem::path& filename) {
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open mol2 template file: " + filename.string());
    }
    std::ostringstream buffer;
    buffer << input.rdbuf();
    auto molecule = load_mol2_text(buffer.str());
    if (molecule.residues.empty()) {
        throw std::invalid_argument("mol2 template file contains no residues: " + filename.string());
    }

    std::vector<ResidueType> residue_types;
    residue_types.reserve(molecule.residues.size());
    for (const auto& residue : molecule.residues) {
        residue_types.push_back(residue_type_from_molecule_residue(molecule, residue));
    }

    for (auto& residue_type : residue_types) {
        put_template(std::move(residue_type));
    }
    std::unique_lock lock(registry_mutex());
    molecule_templates().insert_or_assign(molecule.residues.front().name, std::move(molecule));
}

void register_template_virtual_atom2(const std::string& template_name, const std::string& virtual_atom,
                                     const std::string& atom0, const std::string& atom1, const std::string& atom2,
                                     double k1, double k2) {
    std::unique_lock lock(registry_mutex());
    auto it = molecule_templates().find(template_name);
    if (it == molecule_templates().end()) {
        throw std::out_of_range("molecule template not found: " + template_name);
    }
    auto& molecule = it->second;
    if (molecule.residues.empty()) {
        throw std::invalid_argument("molecule template contains no residues: " + template_name);
    }
    const auto& residue = molecule.residues.front();
    const auto find_atom_by_name = [&](const std::string& name) {
        for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
            const AtomId atom_id = residue.atom_begin + local;
            if (molecule.atoms[atom_id].name == name) {
                return atom_id;
            }
        }
        throw std::out_of_range("atom name not found in molecule template " + template_name + ": " + name);
    };
    molecule.add_virtual_atom2(find_atom_by_name(virtual_atom), find_atom_by_name(atom0),
                               find_atom_by_name(atom1), find_atom_by_name(atom2), k1, k2);
}

void configure_residue_template_head(const std::string& template_name, const std::string& atom,
                                     double length, const std::string& next) {
    std::unique_lock lock(registry_mutex());
    auto it = templates().find(template_name);
    if (it == templates().end()) {
        throw std::out_of_range("residue template not found: " + template_name);
    }
    it->second.set_head(atom, length, next);
}

void configure_residue_template_tail(const std::string& template_name, const std::string& atom,
                                     double length, const std::string& next) {
    std::unique_lock lock(registry_mutex());
    auto it = templates().find(template_name);
    if (it == templates().end()) {
        throw std::out_of_range("residue template not found: " + template_name);
    }
    it->second.set_tail(atom, length, next);
}

void configure_residue_template_connect_atom(const std::string& template_name, const std::string& key,
                                             const std::string& atom) {
    std::unique_lock lock(registry_mutex());
    auto it = templates().find(template_name);
    if (it == templates().end()) {
        throw std::out_of_range("residue template not found: " + template_name);
    }
    it->second.set_connect_atom(key, atom);
}

void register_residue_template_alias(const std::string& alias_name, const std::string& template_name) {
    std::unique_lock lock(registry_mutex());
    const auto it = templates().find(template_name);
    if (it == templates().end()) {
        throw std::out_of_range("residue template not found: " + template_name);
    }
    ResidueType alias(alias_name);
    for (const auto& atom : it->second.atoms()) {
        alias.add_atom(atom.name, atom.type, atom.x, atom.y, atom.z, atom.charge, atom.mass);
    }
    for (const auto& bond : it->second.bonds()) {
        alias.add_connectivity(it->second.atoms()[bond.atom1].name, it->second.atoms()[bond.atom2].name);
    }
    if (!it->second.head().empty()) {
        alias.set_head(it->second.head(), it->second.head_length(), it->second.head_next());
    }
    if (!it->second.tail().empty()) {
        alias.set_tail(it->second.tail(), it->second.tail_length(), it->second.tail_next());
    }
    for (const auto& [key, atom] : it->second.connect_atoms()) {
        alias.set_connect_atom(key, atom);
    }
    templates().insert_or_assign(alias_name, std::move(alias));

    const auto molecule_it = molecule_templates().find(template_name);
    if (molecule_it != molecule_templates().end()) {
        Molecule alias_molecule = molecule_it->second;
        alias_molecule.name = alias_name;
        for (auto& residue : alias_molecule.residues) {
            residue.name = alias_name;
            residue.type_name = alias_name;
            residue.original_name = alias_name;
        }
        molecule_templates().insert_or_assign(alias_name, std::move(alias_molecule));
    }
}

bool has_template(const std::string& name) {
    std::shared_lock lock(registry_mutex());
    return templates().find(name) != templates().end();
}

std::size_t template_atom_count(const std::string& name) {
    std::shared_lock lock(registry_mutex());
    const auto it = templates().find(name);
    if (it == templates().end()) {
        throw std::out_of_range("residue template not found: " + name);
    }
    return it->second.atom_count();
}

std::vector<std::string> registered_template_names() {
    std::shared_lock lock(registry_mutex());
    std::vector<std::string> names;
    names.reserve(templates().size());
    for (const auto& [name, _template] : templates()) {
        names.push_back(name);
    }
    std::sort(names.begin(), names.end());
    return names;
}

const ResidueType& get_residue_template(const std::string& name) {
    std::shared_lock lock(registry_mutex());
    const auto it = templates().find(name);
    if (it == templates().end()) {
        throw std::out_of_range("residue template not found: " + name);
    }
    return it->second;
}

Molecule get_template_molecule(const std::string& name) {
    {
        std::shared_lock lock(registry_mutex());
        const auto molecule_it = molecule_templates().find(name);
        if (molecule_it != molecule_templates().end()) {
            return molecule_it->second;
        }
    }
    const auto& residue_type = get_residue_template(name);
    Molecule molecule(name);
    molecule.append_residue_from_type(residue_type, 0.0, 0.0, 0.0);
    if (!molecule.atoms.empty()) {
        molecule.set_box_padding(0.0, false);
    }
    return molecule;
}

}  // namespace xpongecpp
