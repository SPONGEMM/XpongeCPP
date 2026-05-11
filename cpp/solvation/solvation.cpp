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

ResidueType first_residue_as_type(const Molecule& molecule) {
    if (molecule.residues.empty()) {
        throw std::invalid_argument("solvent molecule must contain one residue");
    }
    const auto& residue = molecule.residues.front();
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
        const auto& atom1 = molecule.atoms[bond.atom1];
        const auto& atom2 = molecule.atoms[bond.atom2];
        residue_type.add_connectivity(atom1.name, atom2.name);
    }
    return residue_type;
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
    const auto solvent_type = first_residue_as_type(solvent);
    const std::size_t atoms_per_solvent = solvent_type.atom_count();

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
    molecule.residues.reserve(molecule.residues.size() + placements.size());
    for (const auto& placement : placements) {
        molecule.append_residue_from_type(solvent_type, placement[0], placement[1], placement[2]);
    }
}

void add_ions(Molecule& molecule, const std::unordered_map<std::string, std::int64_t>& counts, std::uint64_t seed) {
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

    std::vector<ResidueId> water_residues;
    water_residues.reserve(molecule.residues.size());
    for (ResidueId residue_id = 0; residue_id < molecule.residues.size(); ++residue_id) {
        const auto& residue = molecule.residues[residue_id];
        if (residue.name == "WAT") {
            water_residues.push_back(residue_id);
        }
    }
    if (static_cast<std::int64_t>(water_residues.size()) < requested) {
        throw std::invalid_argument("not enough WAT residues to replace with ions");
    }

    std::mt19937_64 rng(seed == 0 ? 5489ULL : seed);
    std::shuffle(water_residues.begin(), water_residues.end(), rng);

    const auto order = deterministic_ion_order(counts);
    std::unordered_map<ResidueId, std::string> replacement_by_residue;
    replacement_by_residue.reserve(static_cast<std::size_t>(requested));
    std::size_t cursor = 0;
    for (const auto& ion_name : order) {
        const auto count = counts.at(ion_name);
        for (std::int64_t i = 0; i < count; ++i) {
            replacement_by_residue.emplace(water_residues[cursor++], ion_name);
        }
    }

    Molecule rebuilt(molecule.name);
    rebuilt.box_length = molecule.box_length;
    rebuilt.box_angle = molecule.box_angle;
    rebuilt.has_box = molecule.has_box;
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

    const auto append_ion_from_water = [&](const Residue& water, const std::string& ion_name) {
        const auto& ion_type = get_residue_template(ion_name);
        const auto& oxygen = molecule.atoms[water.atom_begin];
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
        ion.x = oxygen.x;
        ion.y = oxygen.y;
        ion.z = oxygen.z;
        ion.charge = ion_atom.charge;
        ion.mass = ion_atom.mass;
        rebuilt.atoms.push_back(std::move(ion));
    };

    for (ResidueId residue_id = 0; residue_id < molecule.residues.size(); ++residue_id) {
        const auto& residue = molecule.residues[residue_id];
        if (residue.name != "WAT") {
            append_residue_copy(residue);
        }
    }
    for (const auto& ion_name : order) {
        for (ResidueId residue_id = 0; residue_id < molecule.residues.size(); ++residue_id) {
            const auto it = replacement_by_residue.find(residue_id);
            if (it != replacement_by_residue.end() && it->second == ion_name) {
                append_ion_from_water(molecule.residues[residue_id], ion_name);
            }
        }
    }
    for (ResidueId residue_id = 0; residue_id < molecule.residues.size(); ++residue_id) {
        const auto& residue = molecule.residues[residue_id];
        if (residue.name == "WAT" && replacement_by_residue.find(residue_id) == replacement_by_residue.end()) {
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

    molecule = std::move(rebuilt);
    if (!molecule.validate()) {
        throw std::runtime_error("internal error: invalid molecule after ion replacement");
    }
}

}  // namespace xpongecpp
