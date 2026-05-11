#include "core.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <utility>

namespace xpongecpp {
namespace {

std::array<double, 3> molecule_min(const Molecule& molecule) {
    std::array<double, 3> minv{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
    };
    for (const auto& atom : molecule.atoms) {
        minv[0] = std::min(minv[0], atom.x);
        minv[1] = std::min(minv[1], atom.y);
        minv[2] = std::min(minv[2], atom.z);
    }
    return minv;
}

std::array<double, 3> molecule_max(const Molecule& molecule) {
    std::array<double, 3> maxv{
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
    };
    for (const auto& atom : molecule.atoms) {
        maxv[0] = std::max(maxv[0], atom.x);
        maxv[1] = std::max(maxv[1], atom.y);
        maxv[2] = std::max(maxv[2], atom.z);
    }
    return maxv;
}

std::vector<std::string> deterministic_ion_order(const std::unordered_map<std::string, std::int64_t>& counts) {
    std::vector<std::string> order;
    for (const auto* name : {"NA", "CL"}) {
        const auto it = counts.find(name);
        if (it != counts.end() && it->second > 0) {
            order.emplace_back(name);
        }
    }
    for (const auto& [name, count] : counts) {
        if (name != "NA" && name != "CL" && count > 0) {
            order.push_back(name);
        }
    }
    return order;
}

void append_solvent_unit(Molecule& molecule, const Molecule& solvent, const std::array<double, 3>& placement) {
    const AtomId atom_offset = static_cast<AtomId>(molecule.atoms.size());
    const ResidueId residue_offset = static_cast<ResidueId>(molecule.residues.size());

    molecule.residues.reserve(molecule.residues.size() + solvent.residues.size());
    molecule.atoms.reserve(molecule.atoms.size() + solvent.atoms.size());
    molecule.explicit_bonds.reserve(molecule.explicit_bonds.size() + solvent.explicit_bonds.size());
    molecule.residue_links.reserve(molecule.residue_links.size() + solvent.residue_links.size());
    molecule.virtual_atoms.reserve(molecule.virtual_atoms.size() + solvent.virtual_atoms.size());
    molecule.harmonic_impropers.reserve(molecule.harmonic_impropers.size() + solvent.harmonic_impropers.size());
    molecule.cmap_types.reserve(molecule.cmap_types.size() + solvent.cmap_types.size());
    molecule.cmaps.reserve(molecule.cmaps.size() + solvent.cmaps.size());
    molecule.nb14_extras.reserve(molecule.nb14_extras.size() + solvent.nb14_extras.size());
    molecule.urey_bradleys.reserve(molecule.urey_bradleys.size() + solvent.urey_bradleys.size());
    molecule.ryckaert_bellemans.reserve(molecule.ryckaert_bellemans.size() + solvent.ryckaert_bellemans.size());
    molecule.soft_bonds.reserve(molecule.soft_bonds.size() + solvent.soft_bonds.size());
    molecule.listed_force_definitions.reserve(molecule.listed_force_definitions.size() +
                                              solvent.listed_force_definitions.size());

    for (const auto& residue : solvent.residues) {
        Residue copied = residue;
        copied.atom_begin += atom_offset;
        molecule.residues.push_back(std::move(copied));
    }
    for (const auto& source_atom : solvent.atoms) {
        Atom atom = source_atom;
        atom.residue += residue_offset;
        atom.x += placement[0];
        atom.y += placement[1];
        atom.z += placement[2];
        molecule.atoms.push_back(std::move(atom));
    }
    for (const auto& bond : solvent.explicit_bonds) {
        molecule.explicit_bonds.push_back({bond.atom1 + atom_offset, bond.atom2 + atom_offset});
    }
    for (const auto& link : solvent.residue_links) {
        molecule.residue_links.push_back({link.atom1 + atom_offset, link.atom2 + atom_offset});
    }
    for (const auto& vatom : solvent.virtual_atoms) {
        molecule.virtual_atoms.push_back({vatom.virtual_atom + atom_offset, vatom.atom0 + atom_offset,
                                          vatom.atom1 + atom_offset, vatom.atom2 + atom_offset,
                                          vatom.k1, vatom.k2});
    }
    for (const auto& improper : solvent.harmonic_impropers) {
        molecule.harmonic_impropers.push_back({improper.atom0 + atom_offset, improper.atom1 + atom_offset,
                                               improper.atom2 + atom_offset, improper.atom3 + atom_offset,
                                               improper.k, improper.phi0});
    }
    const std::uint32_t cmap_type_offset = static_cast<std::uint32_t>(molecule.cmap_types.size());
    for (const auto& type : solvent.cmap_types) {
        molecule.cmap_types.push_back(type);
    }
    for (const auto& cmap : solvent.cmaps) {
        molecule.cmaps.push_back({cmap.atom0 + atom_offset, cmap.atom1 + atom_offset, cmap.atom2 + atom_offset,
                                  cmap.atom3 + atom_offset, cmap.atom4 + atom_offset, cmap.type + cmap_type_offset});
    }
    for (const auto& nb14 : solvent.nb14_extras) {
        molecule.nb14_extras.push_back({nb14.atom1 + atom_offset, nb14.atom2 + atom_offset,
                                        nb14.a, nb14.b, nb14.kee});
    }
    for (const auto& angle : solvent.urey_bradleys) {
        molecule.urey_bradleys.push_back({angle.atom0 + atom_offset, angle.atom1 + atom_offset,
                                          angle.atom2 + atom_offset, angle.k, angle.b, angle.k_ub, angle.r13});
    }
    for (const auto& dihedral : solvent.ryckaert_bellemans) {
        molecule.ryckaert_bellemans.push_back({dihedral.atom0 + atom_offset, dihedral.atom1 + atom_offset,
                                               dihedral.atom2 + atom_offset, dihedral.atom3 + atom_offset,
                                               dihedral.c0, dihedral.c1, dihedral.c2,
                                               dihedral.c3, dihedral.c4, dihedral.c5});
    }
    for (const auto& bond : solvent.soft_bonds) {
        molecule.soft_bonds.push_back({bond.atom1 + atom_offset, bond.atom2 + atom_offset,
                                       bond.k, bond.b, bond.from_a_or_b});
    }
    for (const auto& definition : solvent.listed_force_definitions) {
        molecule.listed_force_definitions.push_back(definition);
    }
}

}  // namespace

void add_solvent_box(Molecule& molecule, const Molecule& solvent, double distance, double tolerance,
                     std::int64_t n_solvent, std::uint64_t seed) {
    if (molecule.atoms.empty() || solvent.atoms.empty() || solvent.residues.empty()) {
        throw std::invalid_argument("solute and solvent must contain atoms");
    }
    if (distance < 0.0 || tolerance <= 0.0) {
        throw std::invalid_argument("distance should be non-negative and tolerance should be positive");
    }
    if (n_solvent < 0) {
        throw std::invalid_argument("n_solvent should be non-negative or omitted");
    }

    const auto minv = molecule_min(molecule);
    const auto maxv = molecule_max(molecule);
    const auto solvent_min = molecule_min(solvent);
    const auto solvent_max = molecule_max(solvent);
    const std::array<double, 3> solvent_shape{
        std::max(2.4, solvent_max[0] - solvent_min[0] + tolerance),
        std::max(2.4, solvent_max[1] - solvent_min[1] + tolerance),
        std::max(2.4, solvent_max[2] - solvent_min[2] + tolerance),
    };
    const std::size_t atoms_per_solvent = solvent.atoms.size();
    const std::size_t residues_per_solvent = solvent.residues.size();

    std::vector<std::array<double, 3>> placements;
    const std::array<std::int64_t, 3> inner_grid{
        static_cast<std::int64_t>(std::floor((maxv[0] - minv[0]) / solvent_shape[0])),
        static_cast<std::int64_t>(std::floor((maxv[1] - minv[1]) / solvent_shape[1])),
        static_cast<std::int64_t>(std::floor((maxv[2] - minv[2]) / solvent_shape[2])),
    };
    std::int64_t inner_added = 0;
    if (inner_grid[0] > 0 && inner_grid[1] > 0 && inner_grid[2] > 0) {
        std::vector<std::uint8_t> grid(static_cast<std::size_t>(inner_grid[0] * inner_grid[1] * inner_grid[2]), 1);
        const auto index_of = [&](std::int64_t i, std::int64_t j, std::int64_t k) {
            return static_cast<std::size_t>((i * inner_grid[1] + j) * inner_grid[2] + k);
        };
        for (const auto& atom : molecule.atoms) {
            const std::array<std::int64_t, 3> index{
                static_cast<std::int64_t>(std::floor((atom.x - minv[0]) / solvent_shape[0])),
                static_cast<std::int64_t>(std::floor((atom.y - minv[1]) / solvent_shape[1])),
                static_cast<std::int64_t>(std::floor((atom.z - minv[2]) / solvent_shape[2])),
            };
            for (std::int64_t di = -1; di <= 1; ++di) {
                for (std::int64_t dj = -1; dj <= 1; ++dj) {
                    for (std::int64_t dk = -1; dk <= 1; ++dk) {
                        const auto i = std::clamp(index[0] + di, std::int64_t{0}, inner_grid[0] - 1);
                        const auto j = std::clamp(index[1] + dj, std::int64_t{0}, inner_grid[1] - 1);
                        const auto k = std::clamp(index[2] + dk, std::int64_t{0}, inner_grid[2] - 1);
                        grid[index_of(i, j, k)] = 0;
                    }
                }
            }
        }
        for (std::int64_t i = 0; i < inner_grid[0]; ++i) {
            for (std::int64_t j = 0; j < inner_grid[1]; ++j) {
                for (std::int64_t k = 0; k < inner_grid[2]; ++k) {
                    if (grid[index_of(i, j, k)] != 0) {
                        placements.push_back({minv[0] + i * solvent_shape[0], minv[1] + j * solvent_shape[1],
                                              minv[2] + k * solvent_shape[2]});
                        ++inner_added;
                    }
                }
            }
        }
    }

    const std::array<double, 3> boxmin{minv[0] - distance, minv[1] - distance, minv[2] - distance};
    const std::array<std::int64_t, 3> outer_grid{
        static_cast<std::int64_t>(std::ceil((maxv[0] + distance - boxmin[0]) / solvent_shape[0])),
        static_cast<std::int64_t>(std::ceil((maxv[1] + distance - boxmin[1]) / solvent_shape[1])),
        static_cast<std::int64_t>(std::ceil((maxv[2] + distance - boxmin[2]) / solvent_shape[2])),
    };
    std::vector<std::array<std::int64_t, 3>> outer_candidates;
    const std::array<std::int64_t, 3> in_min{
        static_cast<std::int64_t>(std::floor((minv[0] - boxmin[0]) / solvent_shape[0])),
        static_cast<std::int64_t>(std::floor((minv[1] - boxmin[1]) / solvent_shape[1])),
        static_cast<std::int64_t>(std::floor((minv[2] - boxmin[2]) / solvent_shape[2])),
    };
    const std::array<std::int64_t, 3> in_max{
        static_cast<std::int64_t>(std::ceil((maxv[0] - boxmin[0]) / solvent_shape[0])),
        static_cast<std::int64_t>(std::ceil((maxv[1] - boxmin[1]) / solvent_shape[1])),
        static_cast<std::int64_t>(std::ceil((maxv[2] - boxmin[2]) / solvent_shape[2])),
    };
    for (std::int64_t i = 0; i < outer_grid[0]; ++i) {
        for (std::int64_t j = 0; j < outer_grid[1]; ++j) {
            for (std::int64_t k = 0; k < outer_grid[2]; ++k) {
                const bool inside = i >= in_min[0] && i < in_max[0] && j >= in_min[1] && j < in_max[1] &&
                                    k >= in_min[2] && k < in_max[2];
                if (!inside) {
                    outer_candidates.push_back({i, j, k});
                }
            }
        }
    }
    const std::int64_t max_solvents = inner_added + static_cast<std::int64_t>(outer_candidates.size());
    if (n_solvent > 0 && n_solvent < inner_added) {
        n_solvent = 0;
    }
    if (n_solvent > max_solvents) {
        n_solvent = 0;
    }
    if (n_solvent > 0) {
        const std::int64_t outer_needed = n_solvent - inner_added;
        if (outer_needed < static_cast<std::int64_t>(outer_candidates.size())) {
            std::mt19937_64 rng(seed == 0 ? 5489ULL : seed);
            std::shuffle(outer_candidates.begin(), outer_candidates.end(), rng);
            outer_candidates.resize(static_cast<std::size_t>(outer_needed));
            std::sort(outer_candidates.begin(), outer_candidates.end());
        }
    }
    for (const auto& index : outer_candidates) {
        placements.push_back({boxmin[0] + index[0] * solvent_shape[0], boxmin[1] + index[1] * solvent_shape[1],
                              boxmin[2] + index[2] * solvent_shape[2]});
    }

    molecule.atoms.reserve(molecule.atoms.size() + placements.size() * atoms_per_solvent);
    molecule.residues.reserve(molecule.residues.size() + placements.size() * residues_per_solvent);
    for (const auto& placement : placements) {
        append_solvent_unit(molecule, solvent, placement);
    }
}

void add_ions(Molecule& molecule, const std::unordered_map<std::string, std::int64_t>& counts, std::uint64_t seed,
              const std::string& solvent_residue) {
    std::int64_t requested = 0;
    for (const auto& [name, count] : counts) {
        if (count < 0) {
            throw std::invalid_argument("ion count should be non-negative: " + name);
        }
        requested += count;
    }
    if (requested == 0) {
        return;
    }

    if (solvent_residue.empty()) {
        throw std::invalid_argument("solvent residue name should not be empty");
    }

    std::vector<ResidueId> solvent_residues;
    solvent_residues.reserve(molecule.residues.size());
    for (ResidueId residue_id = 0; residue_id < molecule.residues.size(); ++residue_id) {
        const auto& residue = molecule.residues[residue_id];
        if (residue.name == solvent_residue) {
            solvent_residues.push_back(residue_id);
        }
    }
    if (static_cast<std::int64_t>(solvent_residues.size()) < requested) {
        throw std::invalid_argument("not enough " + solvent_residue + " residues to replace with ions");
    }

    std::mt19937_64 rng(seed == 0 ? 5489ULL : seed);
    std::shuffle(solvent_residues.begin(), solvent_residues.end(), rng);

    const auto order = deterministic_ion_order(counts);
    std::unordered_map<ResidueId, std::string> replacement_by_residue;
    replacement_by_residue.reserve(static_cast<std::size_t>(requested));
    std::size_t cursor = 0;
    for (const auto& ion_name : order) {
        const auto count = counts.at(ion_name);
        for (std::int64_t i = 0; i < count; ++i) {
            replacement_by_residue.emplace(solvent_residues[cursor++], ion_name);
        }
    }

    Molecule rebuilt(molecule.name);
    rebuilt.box_length = molecule.box_length;
    rebuilt.box_angle = molecule.box_angle;
    rebuilt.has_box = molecule.has_box;
    rebuilt.has_gb_parameters = molecule.has_gb_parameters;
    rebuilt.write_min_bonded_parameters = molecule.write_min_bonded_parameters;
    rebuilt.write_subsys_division = molecule.write_subsys_division;
    rebuilt.atoms.reserve(molecule.atoms.size());
    rebuilt.residues.reserve(molecule.residues.size());

    constexpr AtomId invalid_atom_id = std::numeric_limits<AtomId>::max();
    std::vector<AtomId> old_to_new_atom(molecule.atoms.size(), invalid_atom_id);

    const auto append_atom_copy = [&](const Atom& source, AtomId old_atom_id, ResidueId residue_id) {
        Atom atom = source;
        atom.residue = residue_id;
        const AtomId new_atom_id = static_cast<AtomId>(rebuilt.atoms.size());
        rebuilt.atoms.push_back(std::move(atom));
        old_to_new_atom[old_atom_id] = new_atom_id;
    };

    const auto append_residue_copy = [&](const Residue& residue) {
        const ResidueId new_residue_id = static_cast<ResidueId>(rebuilt.residues.size());
        Residue new_residue;
        new_residue.atom_begin = static_cast<AtomId>(rebuilt.atoms.size());
        new_residue.name = residue.name;
        new_residue.type_name = residue.type_name;
        new_residue.atom_count = residue.atom_count;
        rebuilt.residues.push_back(new_residue);
        for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
            const AtomId old_atom_id = residue.atom_begin + local;
            append_atom_copy(molecule.atoms[old_atom_id], old_atom_id, new_residue_id);
        }
    };

    const auto append_ion_from_solvent = [&](const Residue& solvent_residue_to_replace, const std::string& ion_name) {
        const auto& ion_type = get_residue_template(ion_name);
        const auto& anchor = molecule.atoms[solvent_residue_to_replace.atom_begin];
        const ResidueId new_residue_id = static_cast<ResidueId>(rebuilt.residues.size());
        Residue new_residue;
        new_residue.name = ion_type.name();
        new_residue.type_name = ion_type.name();
        new_residue.atom_begin = static_cast<AtomId>(rebuilt.atoms.size());
        new_residue.atom_count = 1;
        rebuilt.residues.push_back(new_residue);

        Atom ion;
        const auto& ion_atom = ion_type.atoms().front();
        ion.name = ion_atom.name;
        ion.type = ion_atom.type;
        ion.element = ion_atom.element;
        ion.residue = new_residue_id;
        ion.x = anchor.x;
        ion.y = anchor.y;
        ion.z = anchor.z;
        ion.charge = ion_atom.charge;
        ion.mass = ion_atom.mass;
        rebuilt.atoms.push_back(std::move(ion));
    };

    for (ResidueId residue_id = 0; residue_id < molecule.residues.size(); ++residue_id) {
        const auto& residue = molecule.residues[residue_id];
        if (residue.name != solvent_residue) {
            append_residue_copy(residue);
        }
    }
    for (const auto& ion_name : order) {
        for (ResidueId residue_id = 0; residue_id < molecule.residues.size(); ++residue_id) {
            const auto it = replacement_by_residue.find(residue_id);
            if (it != replacement_by_residue.end() && it->second == ion_name) {
                append_ion_from_solvent(molecule.residues[residue_id], ion_name);
            }
        }
    }
    for (ResidueId residue_id = 0; residue_id < molecule.residues.size(); ++residue_id) {
        const auto& residue = molecule.residues[residue_id];
        if (residue.name == solvent_residue &&
            replacement_by_residue.find(residue_id) == replacement_by_residue.end()) {
            append_residue_copy(residue);
        }
    }

    rebuilt.explicit_bonds.reserve(molecule.explicit_bonds.size());
    for (const auto& bond : molecule.explicit_bonds) {
        if (bond.atom1 >= old_to_new_atom.size() || bond.atom2 >= old_to_new_atom.size()) {
            continue;
        }
        const AtomId atom1 = old_to_new_atom[bond.atom1];
        const AtomId atom2 = old_to_new_atom[bond.atom2];
        if (atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        rebuilt.explicit_bonds.push_back({atom1, atom2});
    }

    rebuilt.residue_links.reserve(molecule.residue_links.size());
    for (const auto& link : molecule.residue_links) {
        if (link.atom1 >= old_to_new_atom.size() || link.atom2 >= old_to_new_atom.size()) {
            continue;
        }
        const AtomId atom1 = old_to_new_atom[link.atom1];
        const AtomId atom2 = old_to_new_atom[link.atom2];
        if (atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        rebuilt.residue_links.push_back({atom1, atom2});
    }
    const auto mapped_atom = [&](AtomId old_atom_id) {
        if (old_atom_id >= old_to_new_atom.size()) {
            return invalid_atom_id;
        }
        return old_to_new_atom[old_atom_id];
    };
    for (const auto& vatom : molecule.virtual_atoms) {
        const AtomId virtual_atom = mapped_atom(vatom.virtual_atom);
        const AtomId atom0 = mapped_atom(vatom.atom0);
        const AtomId atom1 = mapped_atom(vatom.atom1);
        const AtomId atom2 = mapped_atom(vatom.atom2);
        if (virtual_atom == invalid_atom_id || atom0 == invalid_atom_id ||
            atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        rebuilt.virtual_atoms.push_back({virtual_atom, atom0, atom1, atom2, vatom.k1, vatom.k2});
    }
    for (const auto& improper : molecule.harmonic_impropers) {
        const AtomId atom0 = mapped_atom(improper.atom0);
        const AtomId atom1 = mapped_atom(improper.atom1);
        const AtomId atom2 = mapped_atom(improper.atom2);
        const AtomId atom3 = mapped_atom(improper.atom3);
        if (atom0 == invalid_atom_id || atom1 == invalid_atom_id ||
            atom2 == invalid_atom_id || atom3 == invalid_atom_id) {
            continue;
        }
        rebuilt.harmonic_impropers.push_back({atom0, atom1, atom2, atom3, improper.k, improper.phi0});
    }
    rebuilt.cmap_types = molecule.cmap_types;
    for (const auto& cmap : molecule.cmaps) {
        const AtomId atom0 = mapped_atom(cmap.atom0);
        const AtomId atom1 = mapped_atom(cmap.atom1);
        const AtomId atom2 = mapped_atom(cmap.atom2);
        const AtomId atom3 = mapped_atom(cmap.atom3);
        const AtomId atom4 = mapped_atom(cmap.atom4);
        if (atom0 == invalid_atom_id || atom1 == invalid_atom_id || atom2 == invalid_atom_id ||
            atom3 == invalid_atom_id || atom4 == invalid_atom_id) {
            continue;
        }
        rebuilt.cmaps.push_back({atom0, atom1, atom2, atom3, atom4, cmap.type});
    }
    for (const auto& nb14 : molecule.nb14_extras) {
        const AtomId atom1 = mapped_atom(nb14.atom1);
        const AtomId atom2 = mapped_atom(nb14.atom2);
        if (atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        rebuilt.nb14_extras.push_back({atom1, atom2, nb14.a, nb14.b, nb14.kee});
    }
    for (const auto& angle : molecule.urey_bradleys) {
        const AtomId atom0 = mapped_atom(angle.atom0);
        const AtomId atom1 = mapped_atom(angle.atom1);
        const AtomId atom2 = mapped_atom(angle.atom2);
        if (atom0 == invalid_atom_id || atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        rebuilt.urey_bradleys.push_back({atom0, atom1, atom2, angle.k, angle.b, angle.k_ub, angle.r13});
    }
    for (const auto& dihedral : molecule.ryckaert_bellemans) {
        const AtomId atom0 = mapped_atom(dihedral.atom0);
        const AtomId atom1 = mapped_atom(dihedral.atom1);
        const AtomId atom2 = mapped_atom(dihedral.atom2);
        const AtomId atom3 = mapped_atom(dihedral.atom3);
        if (atom0 == invalid_atom_id || atom1 == invalid_atom_id ||
            atom2 == invalid_atom_id || atom3 == invalid_atom_id) {
            continue;
        }
        rebuilt.ryckaert_bellemans.push_back({atom0, atom1, atom2, atom3, dihedral.c0, dihedral.c1,
                                              dihedral.c2, dihedral.c3, dihedral.c4, dihedral.c5});
    }
    for (const auto& bond : molecule.soft_bonds) {
        const AtomId atom1 = mapped_atom(bond.atom1);
        const AtomId atom2 = mapped_atom(bond.atom2);
        if (atom1 == invalid_atom_id || atom2 == invalid_atom_id) {
            continue;
        }
        rebuilt.soft_bonds.push_back({atom1, atom2, bond.k, bond.b, bond.from_a_or_b});
    }
    rebuilt.listed_force_definitions = molecule.listed_force_definitions;

    molecule = std::move(rebuilt);
    if (!molecule.validate()) {
        throw std::runtime_error("internal error: invalid molecule after ion replacement");
    }
}

}  // namespace xpongecpp
