#include "core.hpp"

#include <algorithm>
#include <limits>
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
    const auto minv = molecule_min(molecule);
    for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
        const auto& atom = molecule.atoms[residue.atom_begin + local];
        residue_type.add_atom(atom.name, atom.type, atom.x - minv[0], atom.y - minv[1], atom.z - minv[2],
                              atom.charge, atom.mass);
    }
    if (residue_type.atom_count() >= 3) {
        residue_type.add_connectivity("O", "H1");
        residue_type.add_connectivity("O", "H2");
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

void append_atom_copy(Molecule& target, const Atom& source, ResidueId residue_id) {
    Atom atom = source;
    atom.residue = residue_id;
    target.atoms.push_back(std::move(atom));
}

}  // namespace

void add_solvent_box(Molecule& molecule, const Molecule& solvent, double distance, double tolerance,
                     std::int64_t n_solvent) {
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
    if (n_solvent == 0) {
        const double lx = maxv[0] - minv[0] + 2.0 * distance;
        const double ly = maxv[1] - minv[1] + 2.0 * distance;
        const double lz = maxv[2] - minv[2] + 2.0 * distance;
        n_solvent = static_cast<std::int64_t>(std::max(1.0, (lx * ly * lz) / 120.0));
    }

    const auto solvent_type = first_residue_as_type(solvent);
    const std::size_t atoms_per_solvent = solvent_type.atom_count();
    molecule.atoms.reserve(molecule.atoms.size() + static_cast<std::size_t>(n_solvent) * atoms_per_solvent);
    molecule.residues.reserve(molecule.residues.size() + static_cast<std::size_t>(n_solvent));

    const std::int64_t nx = std::max<std::int64_t>(1, static_cast<std::int64_t>((maxv[0] - minv[0] + 2 * distance) / solvent_shape[0]));
    const std::int64_t ny = std::max<std::int64_t>(1, static_cast<std::int64_t>((maxv[1] - minv[1] + 2 * distance) / solvent_shape[1]));
    for (std::int64_t i = 0; i < n_solvent; ++i) {
        const double x = minv[0] - distance + static_cast<double>(i % nx) * solvent_shape[0];
        const double y = minv[1] - distance + static_cast<double>((i / nx) % ny) * solvent_shape[1];
        const double z = maxv[2] + distance + static_cast<double>(i / (nx * ny)) * solvent_shape[2];
        molecule.append_residue_from_type(solvent_type, x, y, z);
    }
    molecule.set_box_padding(distance, false);
}

void add_ions(Molecule& molecule, const std::unordered_map<std::string, std::int64_t>& counts) {
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

    std::int64_t water_count = 0;
    for (const auto& residue : molecule.residues) {
        if (residue.name == "WAT") {
            ++water_count;
        }
    }
    if (water_count < requested) {
        throw std::invalid_argument("not enough WAT residues to replace with ions");
    }

    Molecule rebuilt(molecule.name);
    rebuilt.box_length = molecule.box_length;
    rebuilt.box_angle = molecule.box_angle;
    rebuilt.atoms.reserve(molecule.atoms.size());
    rebuilt.residues.reserve(molecule.residues.size());

    auto order = deterministic_ion_order(counts);
    std::size_t ion_kind = 0;
    std::int64_t left_for_kind = order.empty() ? 0 : counts.at(order[ion_kind]);

    for (const auto& residue : molecule.residues) {
        const bool replace_water = residue.name == "WAT" && ion_kind < order.size() && left_for_kind > 0;
        const ResidueId new_residue_id = static_cast<ResidueId>(rebuilt.residues.size());
        Residue new_residue;
        new_residue.atom_begin = static_cast<AtomId>(rebuilt.atoms.size());
        if (replace_water) {
            const auto& ion_type = get_residue_template(order[ion_kind]);
            const auto& oxygen = molecule.atoms[residue.atom_begin];
            new_residue.name = ion_type.name();
            new_residue.type_name = ion_type.name();
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
            --left_for_kind;
            while (left_for_kind == 0 && ++ion_kind < order.size()) {
                left_for_kind = counts.at(order[ion_kind]);
            }
        } else {
            new_residue.name = residue.name;
            new_residue.type_name = residue.type_name;
            new_residue.atom_count = residue.atom_count;
            rebuilt.residues.push_back(new_residue);
            for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
                append_atom_copy(rebuilt, molecule.atoms[residue.atom_begin + local], new_residue_id);
            }
        }
    }

    molecule = std::move(rebuilt);
    if (!molecule.validate()) {
        throw std::runtime_error("internal error: invalid molecule after ion replacement");
    }
}

}  // namespace xpongecpp
