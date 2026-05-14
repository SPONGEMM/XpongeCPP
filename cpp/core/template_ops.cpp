#include "core.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace xpongecpp {
namespace {

constexpr AtomId invalid_atom_id = std::numeric_limits<AtomId>::max();

AtomId remap_atom_id(const std::vector<AtomId>& old_to_new_atom, AtomId old_atom_id) {
    if (old_atom_id >= old_to_new_atom.size()) {
        return invalid_atom_id;
    }
    return old_to_new_atom[old_atom_id];
}

void append_internal_structures(Molecule& target, const Molecule& source, AtomId atom_offset) {
    target.explicit_bonds.reserve(target.explicit_bonds.size() + source.explicit_bonds.size());
    for (const auto& bond : source.explicit_bonds) {
        target.explicit_bonds.push_back({bond.atom1 + atom_offset, bond.atom2 + atom_offset});
    }
    target.residue_links.reserve(target.residue_links.size() + source.residue_links.size());
    for (const auto& link : source.residue_links) {
        target.residue_links.push_back({link.atom1 + atom_offset, link.atom2 + atom_offset});
    }
    target.virtual_atoms.reserve(target.virtual_atoms.size() + source.virtual_atoms.size());
    for (const auto& vatom : source.virtual_atoms) {
        target.virtual_atoms.push_back({vatom.virtual_atom + atom_offset, vatom.atom0 + atom_offset,
                                        vatom.atom1 + atom_offset, vatom.atom2 + atom_offset,
                                        vatom.k1, vatom.k2});
    }
    target.harmonic_impropers.reserve(target.harmonic_impropers.size() + source.harmonic_impropers.size());
    for (const auto& improper : source.harmonic_impropers) {
        target.harmonic_impropers.push_back({improper.atom0 + atom_offset, improper.atom1 + atom_offset,
                                             improper.atom2 + atom_offset, improper.atom3 + atom_offset,
                                             improper.k, improper.phi0});
    }
    const std::uint32_t cmap_type_offset = static_cast<std::uint32_t>(target.cmap_types.size());
    target.cmap_types.reserve(target.cmap_types.size() + source.cmap_types.size());
    for (const auto& type : source.cmap_types) {
        target.cmap_types.push_back(type);
    }
    target.cmaps.reserve(target.cmaps.size() + source.cmaps.size());
    for (const auto& cmap : source.cmaps) {
        target.cmaps.push_back({cmap.atom0 + atom_offset, cmap.atom1 + atom_offset, cmap.atom2 + atom_offset,
                                cmap.atom3 + atom_offset, cmap.atom4 + atom_offset, cmap.type + cmap_type_offset});
    }
    target.nb14_extras.reserve(target.nb14_extras.size() + source.nb14_extras.size());
    for (const auto& nb14 : source.nb14_extras) {
        target.nb14_extras.push_back({nb14.atom1 + atom_offset, nb14.atom2 + atom_offset,
                                      nb14.a, nb14.b, nb14.kee});
    }
    target.urey_bradleys.reserve(target.urey_bradleys.size() + source.urey_bradleys.size());
    for (const auto& angle : source.urey_bradleys) {
        target.urey_bradleys.push_back({angle.atom0 + atom_offset, angle.atom1 + atom_offset,
                                        angle.atom2 + atom_offset, angle.k, angle.b, angle.k_ub, angle.r13});
    }
    target.ryckaert_bellemans.reserve(target.ryckaert_bellemans.size() + source.ryckaert_bellemans.size());
    for (const auto& dihedral : source.ryckaert_bellemans) {
        target.ryckaert_bellemans.push_back({dihedral.atom0 + atom_offset, dihedral.atom1 + atom_offset,
                                             dihedral.atom2 + atom_offset, dihedral.atom3 + atom_offset,
                                             dihedral.c0, dihedral.c1, dihedral.c2,
                                             dihedral.c3, dihedral.c4, dihedral.c5});
    }
    target.soft_bonds.reserve(target.soft_bonds.size() + source.soft_bonds.size());
    for (const auto& bond : source.soft_bonds) {
        target.soft_bonds.push_back({bond.atom1 + atom_offset, bond.atom2 + atom_offset,
                                     bond.k, bond.b, bond.from_a_or_b});
    }
    target.listed_force_definitions.reserve(target.listed_force_definitions.size() +
                                            source.listed_force_definitions.size());
    for (const auto& definition : source.listed_force_definitions) {
        target.listed_force_definitions.push_back(definition);
    }
    for (const auto& [name, parameter] : source.sw_parameters) {
        target.sw_parameters[name] = parameter;
    }
    for (const auto& [name, parameter] : source.edip_parameters) {
        target.edip_parameters[name] = parameter;
    }
}

void remap_internal_structures(const Molecule& source, Molecule& target, const std::vector<AtomId>& old_to_new_atom) {
    target.explicit_bonds.reserve(target.explicit_bonds.size() + source.explicit_bonds.size());
    for (const auto& bond : source.explicit_bonds) {
        const AtomId atom1 = remap_atom_id(old_to_new_atom, bond.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, bond.atom2);
        if (atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        target.explicit_bonds.push_back({atom1, atom2});
    }
    target.residue_links.reserve(target.residue_links.size() + source.residue_links.size());
    for (const auto& link : source.residue_links) {
        const AtomId atom1 = remap_atom_id(old_to_new_atom, link.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, link.atom2);
        if (atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        target.residue_links.push_back({atom1, atom2});
    }
    target.virtual_atoms.reserve(target.virtual_atoms.size() + source.virtual_atoms.size());
    for (const auto& vatom : source.virtual_atoms) {
        const AtomId virtual_atom = remap_atom_id(old_to_new_atom, vatom.virtual_atom);
        const AtomId atom0 = remap_atom_id(old_to_new_atom, vatom.atom0);
        const AtomId atom1 = remap_atom_id(old_to_new_atom, vatom.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, vatom.atom2);
        if (virtual_atom == invalid_atom_id || atom0 == invalid_atom_id ||
            atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        target.virtual_atoms.push_back({virtual_atom, atom0, atom1, atom2, vatom.k1, vatom.k2});
    }
    target.harmonic_impropers.reserve(target.harmonic_impropers.size() + source.harmonic_impropers.size());
    for (const auto& improper : source.harmonic_impropers) {
        const AtomId atom0 = remap_atom_id(old_to_new_atom, improper.atom0);
        const AtomId atom1 = remap_atom_id(old_to_new_atom, improper.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, improper.atom2);
        const AtomId atom3 = remap_atom_id(old_to_new_atom, improper.atom3);
        if (atom0 == invalid_atom_id || atom1 == invalid_atom_id ||
            atom2 == invalid_atom_id || atom3 == invalid_atom_id) {
            continue;
        }
        target.harmonic_impropers.push_back({atom0, atom1, atom2, atom3, improper.k, improper.phi0});
    }
    target.cmap_types = source.cmap_types;
    target.cmaps.reserve(target.cmaps.size() + source.cmaps.size());
    for (const auto& cmap : source.cmaps) {
        const AtomId atom0 = remap_atom_id(old_to_new_atom, cmap.atom0);
        const AtomId atom1 = remap_atom_id(old_to_new_atom, cmap.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, cmap.atom2);
        const AtomId atom3 = remap_atom_id(old_to_new_atom, cmap.atom3);
        const AtomId atom4 = remap_atom_id(old_to_new_atom, cmap.atom4);
        if (atom0 == invalid_atom_id || atom1 == invalid_atom_id || atom2 == invalid_atom_id ||
            atom3 == invalid_atom_id || atom4 == invalid_atom_id) {
            continue;
        }
        target.cmaps.push_back({atom0, atom1, atom2, atom3, atom4, cmap.type});
    }
    target.nb14_extras.reserve(target.nb14_extras.size() + source.nb14_extras.size());
    for (const auto& nb14 : source.nb14_extras) {
        const AtomId atom1 = remap_atom_id(old_to_new_atom, nb14.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, nb14.atom2);
        if (atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        target.nb14_extras.push_back({atom1, atom2, nb14.a, nb14.b, nb14.kee});
    }
    target.urey_bradleys.reserve(target.urey_bradleys.size() + source.urey_bradleys.size());
    for (const auto& angle : source.urey_bradleys) {
        const AtomId atom0 = remap_atom_id(old_to_new_atom, angle.atom0);
        const AtomId atom1 = remap_atom_id(old_to_new_atom, angle.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, angle.atom2);
        if (atom0 == invalid_atom_id || atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        target.urey_bradleys.push_back({atom0, atom1, atom2, angle.k, angle.b, angle.k_ub, angle.r13});
    }
    target.ryckaert_bellemans.reserve(target.ryckaert_bellemans.size() + source.ryckaert_bellemans.size());
    for (const auto& dihedral : source.ryckaert_bellemans) {
        const AtomId atom0 = remap_atom_id(old_to_new_atom, dihedral.atom0);
        const AtomId atom1 = remap_atom_id(old_to_new_atom, dihedral.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, dihedral.atom2);
        const AtomId atom3 = remap_atom_id(old_to_new_atom, dihedral.atom3);
        if (atom0 == invalid_atom_id || atom1 == invalid_atom_id ||
            atom2 == invalid_atom_id || atom3 == invalid_atom_id) {
            continue;
        }
        target.ryckaert_bellemans.push_back({atom0, atom1, atom2, atom3, dihedral.c0, dihedral.c1,
                                             dihedral.c2, dihedral.c3, dihedral.c4, dihedral.c5});
    }
    target.soft_bonds.reserve(target.soft_bonds.size() + source.soft_bonds.size());
    for (const auto& bond : source.soft_bonds) {
        const AtomId atom1 = remap_atom_id(old_to_new_atom, bond.atom1);
        const AtomId atom2 = remap_atom_id(old_to_new_atom, bond.atom2);
        if (atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        target.soft_bonds.push_back({atom1, atom2, bond.k, bond.b, bond.from_a_or_b});
    }
}

}  // namespace

void Molecule::replace_residues(const std::unordered_map<ResidueId, Molecule>& replacements,
                                const std::vector<double>& residue_sort_keys, bool sort) {
    for (const auto& [residue_id, replacement] : replacements) {
        if (residue_id >= residues.size()) {
            throw std::out_of_range("ResidueId out of range");
        }
        if (replacement.residue_count() != 1) {
            throw std::invalid_argument("replacement molecules should contain exactly one residue");
        }
        if (!replacement.validate()) {
            throw std::invalid_argument("replacement molecule is invalid");
        }
    }
    if (!residue_sort_keys.empty() && residue_sort_keys.size() != residues.size()) {
        throw std::invalid_argument("residue_sort_keys should match residue count");
    }

    Molecule rebuilt(name);
    rebuilt.box_length = box_length;
    rebuilt.box_angle = box_angle;
    rebuilt.has_box = has_box;
    rebuilt.has_gb_parameters = has_gb_parameters;
    rebuilt.write_min_bonded_parameters = write_min_bonded_parameters;
    rebuilt.write_subsys_division = write_subsys_division;
    rebuilt.write_lj_soft_core = write_lj_soft_core;
    rebuilt.sw_parameters = sw_parameters;
    rebuilt.edip_parameters = edip_parameters;

    std::vector<AtomId> old_to_new_atom(atoms.size(), invalid_atom_id);
    std::vector<ResidueId> residue_order(residues.size());
    for (ResidueId residue_id = 0; residue_id < residues.size(); ++residue_id) {
        residue_order[residue_id] = residue_id;
    }
    if (sort && !residue_sort_keys.empty()) {
        std::stable_sort(residue_order.begin(), residue_order.end(),
                         [&](ResidueId lhs, ResidueId rhs) {
                             return residue_sort_keys[lhs] < residue_sort_keys[rhs];
                         });
    }

    rebuilt.residues.reserve(residues.size());
    for (const ResidueId old_residue_id : residue_order) {
        const auto replacement_it = replacements.find(old_residue_id);
        if (replacement_it == replacements.end()) {
            const auto& source_residue = residues[old_residue_id];
            const ResidueId new_residue_id = static_cast<ResidueId>(rebuilt.residues.size());
            Residue copied_residue = source_residue;
            copied_residue.atom_begin = static_cast<AtomId>(rebuilt.atoms.size());
            rebuilt.residues.push_back(copied_residue);
            for (std::uint32_t local = 0; local < source_residue.atom_count; ++local) {
                const AtomId old_atom_id = source_residue.atom_begin + local;
                Atom atom_copy = atoms[old_atom_id];
                atom_copy.residue = new_residue_id;
                const AtomId new_atom_id = static_cast<AtomId>(rebuilt.atoms.size());
                rebuilt.atoms.push_back(std::move(atom_copy));
                old_to_new_atom[old_atom_id] = new_atom_id;
            }
            continue;
        }

        const auto& source_residue = residues[old_residue_id];
        const auto& replacement = replacement_it->second;
        const auto& replacement_residue = replacement.residues.front();
        const ResidueId new_residue_id = static_cast<ResidueId>(rebuilt.residues.size());
        const AtomId atom_offset = static_cast<AtomId>(rebuilt.atoms.size());
        Residue copied_residue = replacement_residue;
        copied_residue.chain_id = source_residue.chain_id;
        copied_residue.effective_chain_id = source_residue.effective_chain_id;
        copied_residue.segment_id = source_residue.segment_id;
        copied_residue.pdb_resseq = source_residue.pdb_resseq;
        copied_residue.insertion_code = source_residue.insertion_code;
        copied_residue.is_hetero = source_residue.is_hetero;
        copied_residue.atom_begin = atom_offset;
        copied_residue.atom_count = replacement_residue.atom_count;
        rebuilt.residues.push_back(copied_residue);

        double dx = 0.0;
        double dy = 0.0;
        double dz = 0.0;
        if (source_residue.atom_count > 0 && replacement_residue.atom_count > 0) {
            const auto& source_anchor = atoms[source_residue.atom_begin];
            const auto& replacement_anchor = replacement.atoms[replacement_residue.atom_begin];
            dx = source_anchor.x - replacement_anchor.x;
            dy = source_anchor.y - replacement_anchor.y;
            dz = source_anchor.z - replacement_anchor.z;
        }

        for (std::uint32_t local = 0; local < replacement_residue.atom_count; ++local) {
            Atom atom_copy = replacement.atoms[replacement_residue.atom_begin + local];
            atom_copy.residue = new_residue_id;
            atom_copy.x += dx;
            atom_copy.y += dy;
            atom_copy.z += dz;
            rebuilt.atoms.push_back(std::move(atom_copy));
        }
        append_internal_structures(rebuilt, replacement, atom_offset);
    }

    remap_internal_structures(*this, rebuilt, old_to_new_atom);
    *this = std::move(rebuilt);
    if (!validate()) {
        throw std::runtime_error("internal error: invalid molecule after residue replacement");
    }
}

void Molecule::reorder_atoms_by_template(const Molecule& template_molecule) {
    if (residue_count() != template_molecule.residue_count()) {
        throw std::invalid_argument("template molecule should have the same residue count");
    }

    Molecule rebuilt(name);
    rebuilt.box_length = box_length;
    rebuilt.box_angle = box_angle;
    rebuilt.has_box = has_box;
    rebuilt.has_gb_parameters = has_gb_parameters;
    rebuilt.write_min_bonded_parameters = write_min_bonded_parameters;
    rebuilt.write_subsys_division = write_subsys_division;
    rebuilt.write_lj_soft_core = write_lj_soft_core;

    std::vector<AtomId> old_to_new_atom(atoms.size(), invalid_atom_id);
    rebuilt.residues.reserve(residues.size());
    rebuilt.atoms.reserve(atoms.size());

    for (ResidueId residue_id = 0; residue_id < residues.size(); ++residue_id) {
        const auto& source_residue = residues[residue_id];
        const auto& template_residue = template_molecule.residues[residue_id];
        const auto& source_type = !source_residue.type_name.empty() ? source_residue.type_name : source_residue.name;
        const auto& template_type = !template_residue.type_name.empty() ? template_residue.type_name : template_residue.name;
        if (source_type != template_type) {
            throw std::invalid_argument("residue types should match when sorting atoms by template");
        }
        if (source_residue.atom_count != template_residue.atom_count) {
            throw std::invalid_argument("residue atom counts should match when sorting atoms by template");
        }

        std::unordered_map<std::string, AtomId> source_name_to_atom;
        source_name_to_atom.reserve(source_residue.atom_count);
        for (std::uint32_t local = 0; local < source_residue.atom_count; ++local) {
            const AtomId atom_id = source_residue.atom_begin + local;
            const auto inserted = source_name_to_atom.emplace(atoms[atom_id].name, atom_id);
            if (!inserted.second) {
                throw std::invalid_argument("duplicate atom names are not supported in sort_atoms_by");
            }
        }

        Residue copied_residue = source_residue;
        copied_residue.atom_begin = static_cast<AtomId>(rebuilt.atoms.size());
        rebuilt.residues.push_back(copied_residue);
        for (std::uint32_t local = 0; local < template_residue.atom_count; ++local) {
            const auto& template_atom = template_molecule.atoms[template_residue.atom_begin + local];
            const auto found = source_name_to_atom.find(template_atom.name);
            if (found == source_name_to_atom.end()) {
                throw std::invalid_argument("template atom name not found in source residue: " + template_atom.name);
            }
            const AtomId old_atom_id = found->second;
            Atom atom_copy = atoms[old_atom_id];
            atom_copy.residue = residue_id;
            const AtomId new_atom_id = static_cast<AtomId>(rebuilt.atoms.size());
            rebuilt.atoms.push_back(std::move(atom_copy));
            old_to_new_atom[old_atom_id] = new_atom_id;
        }
    }

    remap_internal_structures(*this, rebuilt, old_to_new_atom);
    rebuilt.listed_force_definitions = listed_force_definitions;
    rebuilt.sw_parameters = sw_parameters;
    rebuilt.edip_parameters = edip_parameters;
    *this = std::move(rebuilt);
    if (!validate()) {
        throw std::runtime_error("internal error: invalid molecule after atom reordering");
    }
}

}  // namespace xpongecpp
