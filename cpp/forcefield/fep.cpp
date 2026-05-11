#include "core.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>

namespace xpongecpp {
namespace {

constexpr const char* kZeroLJ = "ZERO_LJ_ATOM";

std::uint64_t pair_key(AtomId atom1, AtomId atom2) {
    const auto lo = std::min(atom1, atom2);
    const auto hi = std::max(atom1, atom2);
    return (static_cast<std::uint64_t>(lo) << 32U) | hi;
}

bool is_zero_lj(const Atom& atom) {
    return atom.type == kZeroLJ;
}

bool all_active(const Molecule& molecule, const std::vector<AtomId>& atoms) {
    return std::all_of(atoms.begin(), atoms.end(), [&](AtomId atom_id) {
        return atom_id < molecule.atoms.size() && !is_zero_lj(molecule.atoms[atom_id]);
    });
}

double distance(const Atom& atom1, const Atom& atom2) {
    const double dx = atom1.x - atom2.x;
    const double dy = atom1.y - atom2.y;
    const double dz = atom1.z - atom2.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

std::pair<double, double> bond_parameters_or_distance(const Atom& atom1, const Atom& atom2) {
    if (const auto term = find_amber_bond_term(atom1.type, atom2.type)) {
        return {term->k, term->length};
    }
    return {300.0, distance(atom1, atom2)};
}

std::string bond_key(const Bond& bond) {
    return std::to_string(std::min(bond.atom1, bond.atom2)) + ":" + std::to_string(std::max(bond.atom1, bond.atom2));
}

std::string angle_key(const Angle& angle) {
    const auto atom1 = std::min(angle.atom1, angle.atom3);
    const auto atom3 = std::max(angle.atom1, angle.atom3);
    return std::to_string(atom1) + ":" + std::to_string(angle.atom2) + ":" + std::to_string(atom3);
}

std::string dihedral_key(const Dihedral& dihedral) {
    std::array<AtomId, 4> forward{dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4};
    std::array<AtomId, 4> reverse{dihedral.atom4, dihedral.atom3, dihedral.atom2, dihedral.atom1};
    const auto& atoms = reverse < forward ? reverse : forward;
    return std::to_string(atoms[0]) + ":" + std::to_string(atoms[1]) + ":" + std::to_string(atoms[2]) + ":" +
           std::to_string(atoms[3]) + ":" + std::to_string(dihedral.periodicity) + ":" +
           std::to_string(dihedral.phase);
}

std::string nb14_key(const NB14& nb14) {
    return std::to_string(std::min(nb14.atom1, nb14.atom2)) + ":" + std::to_string(std::max(nb14.atom1, nb14.atom2));
}

template <typename T, typename KeyFn>
std::unordered_map<std::string, T> active_force_map(const std::vector<T>& forces, const Molecule& molecule,
                                                    KeyFn key_fn) {
    std::unordered_map<std::string, T> out;
    for (const auto& force : forces) {
        if constexpr (std::is_same_v<T, Bond>) {
            if (!all_active(molecule, {force.atom1, force.atom2})) continue;
        } else if constexpr (std::is_same_v<T, Angle>) {
            if (!all_active(molecule, {force.atom1, force.atom2, force.atom3})) continue;
        } else if constexpr (std::is_same_v<T, Dihedral>) {
            if (!all_active(molecule, {force.atom1, force.atom2, force.atom3, force.atom4})) continue;
        } else if constexpr (std::is_same_v<T, NB14>) {
            if (!all_active(molecule, {force.atom1, force.atom2})) continue;
        }
        out.insert_or_assign(key_fn(force), force);
    }
    return out;
}

Topology merge_topology_for_fep(const Molecule& molecule_a, const Molecule& molecule_b,
                                double default_lambda, const std::unordered_map<std::string, double>& specific_lambda) {
    const auto lambda_for = [&](const std::string& key) {
        const auto it = specific_lambda.find(key);
        return it == specific_lambda.end() ? default_lambda : it->second;
    };
    const double bond_lambda = lambda_for("bond_base");
    const double angle_lambda = lambda_for("angle");
    const double dihedral_lambda = lambda_for("dihedral");

    const auto topo_a = build_topology(molecule_a);
    const auto topo_b = build_topology(molecule_b);
    Topology out;

    const auto bonds_a = active_force_map(topo_a.bonds, molecule_a, bond_key);
    const auto bonds_b = active_force_map(topo_b.bonds, molecule_b, bond_key);
    std::unordered_set<std::string> keys;
    for (const auto& [key, force] : bonds_a) keys.insert(key);
    for (const auto& [key, force] : bonds_b) keys.insert(key);
    for (const auto& key : keys) {
        const auto ia = bonds_a.find(key);
        const auto ib = bonds_b.find(key);
        if (ia != bonds_a.end() && ib != bonds_b.end() && std::abs(ia->second.length - ib->second.length) < 1e-5) {
            auto force = ia->second;
            force.k = ia->second.k * (1.0 - bond_lambda) + ib->second.k * bond_lambda;
            out.bonds.push_back(force);
        } else {
            if (ia != bonds_a.end()) {
                auto force = ia->second;
                force.k *= (1.0 - bond_lambda);
                out.bonds.push_back(force);
            }
            if (ib != bonds_b.end()) {
                auto force = ib->second;
                force.k *= bond_lambda;
                out.bonds.push_back(force);
            }
        }
    }

    const auto angles_a = active_force_map(topo_a.angles, molecule_a, angle_key);
    const auto angles_b = active_force_map(topo_b.angles, molecule_b, angle_key);
    keys.clear();
    for (const auto& [key, force] : angles_a) keys.insert(key);
    for (const auto& [key, force] : angles_b) keys.insert(key);
    for (const auto& key : keys) {
        const auto ia = angles_a.find(key);
        const auto ib = angles_b.find(key);
        if (ia != angles_a.end() && ib != angles_b.end() && std::abs(ia->second.theta - ib->second.theta) < 1e-5) {
            auto force = ia->second;
            force.k = ia->second.k * (1.0 - angle_lambda) + ib->second.k * angle_lambda;
            out.angles.push_back(force);
        } else {
            if (ia != angles_a.end()) {
                auto force = ia->second;
                force.k *= (1.0 - angle_lambda);
                out.angles.push_back(force);
            }
            if (ib != angles_b.end()) {
                auto force = ib->second;
                force.k *= angle_lambda;
                out.angles.push_back(force);
            }
        }
    }

    const auto dihedrals_a = active_force_map(topo_a.dihedrals, molecule_a, dihedral_key);
    const auto dihedrals_b = active_force_map(topo_b.dihedrals, molecule_b, dihedral_key);
    keys.clear();
    for (const auto& [key, force] : dihedrals_a) keys.insert(key);
    for (const auto& [key, force] : dihedrals_b) keys.insert(key);
    for (const auto& key : keys) {
        const auto ia = dihedrals_a.find(key);
        const auto ib = dihedrals_b.find(key);
        if (ia != dihedrals_a.end() && ib != dihedrals_b.end()) {
            auto force = ia->second;
            force.k = ia->second.k * (1.0 - dihedral_lambda) + ib->second.k * dihedral_lambda;
            out.dihedrals.push_back(force);
        } else {
            if (ia != dihedrals_a.end()) {
                auto force = ia->second;
                force.k *= (1.0 - dihedral_lambda);
                out.dihedrals.push_back(force);
            }
            if (ib != dihedrals_b.end()) {
                auto force = ib->second;
                force.k *= dihedral_lambda;
                out.dihedrals.push_back(force);
            }
        }
    }

    const auto nb14_a = active_force_map(topo_a.nb14s, molecule_a, nb14_key);
    const auto nb14_b = active_force_map(topo_b.nb14s, molecule_b, nb14_key);
    keys.clear();
    for (const auto& [key, force] : nb14_a) keys.insert(key);
    for (const auto& [key, force] : nb14_b) keys.insert(key);
    for (const auto& key : keys) {
        const auto ia = nb14_a.find(key);
        const auto ib = nb14_b.find(key);
        if (ia != nb14_a.end() && ib != nb14_b.end()) {
            auto force = ia->second;
            force.k_lj = ia->second.k_lj * (1.0 - dihedral_lambda) + ib->second.k_lj * dihedral_lambda;
            force.k_ee = ia->second.k_ee * (1.0 - dihedral_lambda) + ib->second.k_ee * dihedral_lambda;
            out.nb14s.push_back(force);
        } else if (ia != nb14_a.end()) {
            out.nb14s.push_back(ia->second);
        } else if (ib != nb14_b.end()) {
            out.nb14s.push_back(ib->second);
        }
    }

    out.exclusions.resize(molecule_a.atoms.size());
    for (std::size_t i = 0; i < out.exclusions.size(); ++i) {
        std::unordered_set<AtomId> exclusion_set;
        if (i < topo_a.exclusions.size()) {
            exclusion_set.insert(topo_a.exclusions[i].begin(), topo_a.exclusions[i].end());
        }
        if (i < topo_b.exclusions.size()) {
            exclusion_set.insert(topo_b.exclusions[i].begin(), topo_b.exclusions[i].end());
        }
        out.exclusions[i].assign(exclusion_set.begin(), exclusion_set.end());
        std::sort(out.exclusions[i].begin(), out.exclusions[i].end());
    }

    const auto bond_less = [](const Bond& lhs, const Bond& rhs) {
        return std::tie(lhs.atom1, lhs.atom2, lhs.k, lhs.length) < std::tie(rhs.atom1, rhs.atom2, rhs.k, rhs.length);
    };
    std::sort(out.bonds.begin(), out.bonds.end(), bond_less);
    std::sort(out.angles.begin(), out.angles.end(), [](const Angle& lhs, const Angle& rhs) {
        return std::tie(lhs.atom1, lhs.atom2, lhs.atom3, lhs.k, lhs.theta) <
               std::tie(rhs.atom1, rhs.atom2, rhs.atom3, rhs.k, rhs.theta);
    });
    std::sort(out.dihedrals.begin(), out.dihedrals.end(), [](const Dihedral& lhs, const Dihedral& rhs) {
        return std::tie(lhs.atom1, lhs.atom2, lhs.atom3, lhs.atom4, lhs.periodicity, lhs.k, lhs.phase) <
               std::tie(rhs.atom1, rhs.atom2, rhs.atom3, rhs.atom4, rhs.periodicity, rhs.k, rhs.phase);
    });
    std::sort(out.nb14s.begin(), out.nb14s.end(), [](const NB14& lhs, const NB14& rhs) {
        return std::tie(lhs.atom1, lhs.atom2) < std::tie(rhs.atom1, rhs.atom2);
    });
    return out;
}

void add_unique_bond(std::vector<ResidueLink>& bonds, std::unordered_set<std::uint64_t>& seen,
                     AtomId atom1, AtomId atom2) {
    if (atom1 == atom2) {
        return;
    }
    if (seen.insert(pair_key(atom1, atom2)).second) {
        bonds.push_back({std::min(atom1, atom2), std::max(atom1, atom2)});
    }
}

Molecule build_dual_state(const Molecule& molecule, ResidueId residue_index, const Molecule& residue_b_molecule,
                          const std::vector<int>& a_to_b, const std::vector<int>& b_to_a, bool from_a) {
    const auto invalid = std::numeric_limits<AtomId>::max();
    const auto& residue_a = molecule.residue(residue_index);
    const auto& residue_b = residue_b_molecule.residue(0);

    Molecule out(molecule.name + (from_a ? "A" : "B"));
    out.box_length = molecule.box_length;
    out.box_angle = molecule.box_angle;
    out.has_box = molecule.has_box;

    std::vector<AtomId> old_to_new(molecule.atoms.size(), invalid);
    std::vector<AtomId> combined_b_to_new(residue_b.atom_count, invalid);
    std::unordered_set<std::uint64_t> seen_bonds;

    for (ResidueId rid = 0; rid < molecule.residues.size(); ++rid) {
        const auto& source_residue = molecule.residues[rid];
        Residue residue;
        residue.name = source_residue.name;
        residue.type_name = source_residue.type_name;
        residue.atom_begin = static_cast<AtomId>(out.atoms.size());
        residue.atom_count = 0;

        if (rid != residue_index) {
            for (std::uint32_t local = 0; local < source_residue.atom_count; ++local) {
                const AtomId old_atom_id = source_residue.atom_begin + local;
                Atom atom = molecule.atoms[old_atom_id];
                atom.residue = static_cast<ResidueId>(out.residues.size());
                old_to_new[old_atom_id] = static_cast<AtomId>(out.atoms.size());
                out.atoms.push_back(std::move(atom));
                ++residue.atom_count;
            }
            out.residues.push_back(std::move(residue));
            continue;
        }

        residue.name = from_a ? source_residue.name + "_" + residue_b.name : residue_b.name + "_" + source_residue.name;
        residue.type_name = residue.name;
        const auto new_residue_id = static_cast<ResidueId>(out.residues.size());
        for (std::uint32_t local = 0; local < source_residue.atom_count; ++local) {
            const AtomId old_atom_id = source_residue.atom_begin + local;
            Atom atom = molecule.atoms[old_atom_id];
            const int matched_b = a_to_b[local];
            if (!from_a) {
                if (matched_b >= 0) {
                    const auto& b_atom = residue_b_molecule.atoms[residue_b.atom_begin + static_cast<AtomId>(matched_b)];
                    atom.type = b_atom.type;
                    atom.element = b_atom.element;
                    atom.charge = b_atom.charge;
                    atom.mass = b_atom.mass;
                } else {
                    atom.type = kZeroLJ;
                    atom.charge = 0.0;
                    atom.subsys = 1;
                }
            } else if (matched_b < 0) {
                atom.subsys = 1;
            }
            atom.residue = new_residue_id;
            const AtomId new_atom_id = static_cast<AtomId>(out.atoms.size());
            old_to_new[old_atom_id] = new_atom_id;
            out.atoms.push_back(std::move(atom));
            ++residue.atom_count;
            if (matched_b >= 0) {
                combined_b_to_new[static_cast<std::size_t>(matched_b)] = new_atom_id;
            }
        }

        for (std::uint32_t local = 0; local < residue_b.atom_count; ++local) {
            if (b_to_a[local] >= 0) {
                continue;
            }
            Atom atom = residue_b_molecule.atoms[residue_b.atom_begin + local];
            atom.name += "R2";
            atom.subsys = 2;
            if (from_a) {
                atom.type = kZeroLJ;
                atom.charge = 0.0;
            }
            atom.residue = new_residue_id;
            combined_b_to_new[local] = static_cast<AtomId>(out.atoms.size());
            out.atoms.push_back(std::move(atom));
            ++residue.atom_count;
        }
        out.residues.push_back(std::move(residue));
    }

    for (const auto& bond : molecule.explicit_bonds) {
        if (bond.atom1 < old_to_new.size() && bond.atom2 < old_to_new.size() &&
            old_to_new[bond.atom1] != invalid && old_to_new[bond.atom2] != invalid) {
            add_unique_bond(out.explicit_bonds, seen_bonds, old_to_new[bond.atom1], old_to_new[bond.atom2]);
        }
    }
    for (const auto& bond : residue_b_molecule.explicit_bonds) {
        const auto local1 = bond.atom1 - residue_b.atom_begin;
        const auto local2 = bond.atom2 - residue_b.atom_begin;
        if (local1 < combined_b_to_new.size() && local2 < combined_b_to_new.size() &&
            combined_b_to_new[local1] != invalid && combined_b_to_new[local2] != invalid) {
            add_unique_bond(out.explicit_bonds, seen_bonds, combined_b_to_new[local1], combined_b_to_new[local2]);
        }
    }
    for (const auto& link : molecule.residue_links) {
        if (link.atom1 < old_to_new.size() && link.atom2 < old_to_new.size() &&
            old_to_new[link.atom1] != invalid && old_to_new[link.atom2] != invalid) {
            out.residue_links.push_back({old_to_new[link.atom1], old_to_new[link.atom2]});
        }
    }
    if (!out.validate()) {
        throw std::runtime_error("internal error: invalid dual topology molecule");
    }
    return out;
}

}  // namespace

std::pair<Molecule, Molecule> merge_dual_topology(
    const Molecule& molecule, ResidueId residue_index, const Molecule& residue_b_molecule,
    const std::unordered_map<std::uint32_t, std::uint32_t>& match_b_to_a) {
    if (!molecule.validate() || !residue_b_molecule.validate()) {
        throw std::invalid_argument("cannot merge invalid molecule");
    }
    if (residue_index >= molecule.residues.size()) {
        throw std::out_of_range("FEP residue index out of range");
    }
    if (residue_b_molecule.residues.empty()) {
        throw std::invalid_argument("FEP target molecule has no residue");
    }
    const auto& residue_a = molecule.residues[residue_index];
    const auto& residue_b = residue_b_molecule.residues.front();

    std::vector<int> a_to_b(residue_a.atom_count, -1);
    std::vector<int> b_to_a(residue_b.atom_count, -1);
    for (const auto& [b_local, a_local] : match_b_to_a) {
        if (a_local >= residue_a.atom_count || b_local >= residue_b.atom_count) {
            throw std::out_of_range("FEP atom match index out of range");
        }
        if (a_to_b[a_local] >= 0 || b_to_a[b_local] >= 0) {
            throw std::invalid_argument("FEP atom match should be one-to-one");
        }
        a_to_b[a_local] = static_cast<int>(b_local);
        b_to_a[b_local] = static_cast<int>(a_local);
    }

    return {
        build_dual_state(molecule, residue_index, residue_b_molecule, a_to_b, b_to_a, true),
        build_dual_state(molecule, residue_index, residue_b_molecule, a_to_b, b_to_a, false),
    };
}

Molecule merge_force_field(const Molecule& molecule_a, const Molecule& molecule_b, double default_lambda,
                           const std::unordered_map<std::string, double>& specific_lambda) {
    if (!molecule_a.validate() || !molecule_b.validate()) {
        throw std::invalid_argument("cannot merge invalid FEP molecule");
    }
    if (molecule_a.atoms.size() != molecule_b.atoms.size() ||
        molecule_a.residues.size() != molecule_b.residues.size()) {
        throw std::invalid_argument("FEP molecules should have matching atom and residue counts");
    }
    const auto lambda_for = [&](const std::string& key) {
        const auto it = specific_lambda.find(key);
        return it == specific_lambda.end() ? default_lambda : it->second;
    };
    const double charge_lambda = lambda_for("charge");
    const double mass_lambda = lambda_for("mass");

    Molecule out = molecule_a;
    out.name = molecule_a.name + "_merged";
    out.topology_override = std::nullopt;
    for (AtomId atom_id = 0; atom_id < out.atoms.size(); ++atom_id) {
        auto& atom = out.atoms[atom_id];
        const auto& atom_a = molecule_a.atoms[atom_id];
        const auto& atom_b = molecule_b.atoms[atom_id];
        atom.charge = atom_a.charge * (1.0 - charge_lambda) + atom_b.charge * charge_lambda;
        atom.mass = atom_a.mass * (1.0 - mass_lambda) + atom_b.mass * mass_lambda;
        atom.lj_type_b = atom_b.type;
        atom.subsys = atom_a.subsys != 0 ? atom_a.subsys : atom_b.subsys;
    }
    out.enable_lj_soft_core();
    out.topology_override = merge_topology_for_fep(molecule_a, molecule_b, default_lambda, specific_lambda);

    std::unordered_set<std::uint64_t> seen_soft_bonds;
    out.soft_bonds.clear();
    for (const auto& bond : out.explicit_bonds) {
        const bool active_a = !is_zero_lj(molecule_a.atoms[bond.atom1]) && !is_zero_lj(molecule_a.atoms[bond.atom2]);
        const bool active_b = !is_zero_lj(molecule_b.atoms[bond.atom1]) && !is_zero_lj(molecule_b.atoms[bond.atom2]);
        if (active_a == active_b) {
            continue;
        }
        if (!seen_soft_bonds.insert(pair_key(bond.atom1, bond.atom2)).second) {
            continue;
        }
        if (active_a) {
            const auto [k, b] = bond_parameters_or_distance(molecule_a.atoms[bond.atom1], molecule_a.atoms[bond.atom2]);
            out.add_bond_soft(bond.atom1, bond.atom2, k, b, 0);
        } else {
            const auto [k, b] = bond_parameters_or_distance(molecule_b.atoms[bond.atom1], molecule_b.atoms[bond.atom2]);
            out.add_bond_soft(bond.atom1, bond.atom2, k, b, 1);
        }
    }

    return out;
}

}  // namespace xpongecpp
