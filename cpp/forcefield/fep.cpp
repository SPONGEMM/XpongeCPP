#include "core.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <iomanip>
#include <map>
#include <sstream>
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
    return atom.type == kZeroLJ && !atom.fep_preserve_bonded;
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

std::string exact_number(double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("non-finite FEP parameter");
    std::ostringstream out;
    out << std::hexfloat << (value == 0.0 ? 0.0 : value);
    return out.str();
}

std::string bond_key(const Bond& bond) {
    return std::to_string(std::min(bond.atom1, bond.atom2)) + ":" + std::to_string(std::max(bond.atom1, bond.atom2)) + ":" + exact_number(bond.length);
}

std::string angle_key(const Angle& angle) {
    const auto atom1 = std::min(angle.atom1, angle.atom3);
    const auto atom3 = std::max(angle.atom1, angle.atom3);
    return std::to_string(atom1) + ":" + std::to_string(angle.atom2) + ":" + std::to_string(atom3) + ":" + exact_number(angle.theta);
}

std::string dihedral_key(const Dihedral& dihedral) {
    std::array<AtomId, 4> forward{dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4};
    std::array<AtomId, 4> reverse{dihedral.atom4, dihedral.atom3, dihedral.atom2, dihedral.atom1};
    const auto& atoms = reverse < forward ? reverse : forward;
    return std::to_string(atoms[0]) + ":" + std::to_string(atoms[1]) + ":" + std::to_string(atoms[2]) + ":" +
           std::to_string(atoms[3]) + ":" + std::to_string(dihedral.periodicity) + ":" +
           exact_number(dihedral.phase);
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
        // A force is a sum of terms, not a last-value parameter registry.
        // Only identical functions may be coalesced (never arbitrary improper permutations).
        auto [it, inserted] = out.emplace(key_fn(force), force);
        if (!inserted) {
            if constexpr (std::is_same_v<T, NB14>) {
                it->second.k_lj += force.k_lj;
                it->second.k_ee += force.k_ee;
            } else {
                it->second.k += force.k;
            }
        }
    }
    return out;
}

Topology merge_topology_for_fep(const Molecule& molecule_a, const Molecule& molecule_b,
                                double default_lambda, const std::unordered_map<std::string, double>& specific_lambda,
                                std::vector<SoftBond>& soft_bonds) {
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
    std::unordered_set<std::uint64_t> pairs_a, pairs_b;
    for (const auto& [key, force] : bonds_a) pairs_a.insert(pair_key(force.atom1, force.atom2));
    for (const auto& [key, force] : bonds_b) pairs_b.insert(pair_key(force.atom1, force.atom2));
    std::unordered_set<std::string> keys;
    for (const auto& [key, force] : bonds_a) keys.insert(key);
    for (const auto& [key, force] : bonds_b) keys.insert(key);
    for (const auto& key : keys) {
        const auto ia = bonds_a.find(key);
        const auto ib = bonds_b.find(key);
        // A bond that exists in only one physical state is evaluated by
        // SPONGE's soft-bond force. Do not also add a lambda-weighted ordinary
        // harmonic term: that would double count it and leave spurious hard
        // components at the endpoints. Preserve the actual endpoint parameters.
        if (ia != bonds_a.end() && !pairs_b.count(pair_key(ia->second.atom1, ia->second.atom2))) {
            const auto& force = ia->second;
            soft_bonds.push_back({force.atom1, force.atom2, force.k, force.length, 0});
            continue;
        }
        if (ib != bonds_b.end() && !pairs_a.count(pair_key(ib->second.atom1, ib->second.atom2))) {
            const auto& force = ib->second;
            soft_bonds.push_back({force.atom1, force.atom2, force.k, force.length, 1});
            continue;
        }
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
            auto force = ia->second;
            force.k_lj *= (1.0 - dihedral_lambda);
            force.k_ee *= (1.0 - dihedral_lambda);
            out.nb14s.push_back(force);
        } else if (ib != nb14_b.end()) {
            auto force = ib->second;
            force.k_lj *= dihedral_lambda;
            force.k_ee *= dihedral_lambda;
            out.nb14s.push_back(force);
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

std::pair<double, double> fep_lj_ab(const Molecule& molecule, AtomId i, AtomId j) {
    const auto li = resolve_molecule_lj_parameter(molecule, resolve_molecule_lj_type(molecule, molecule.atoms.at(i).type));
    const auto lj = resolve_molecule_lj_parameter(molecule, resolve_molecule_lj_type(molecule, molecule.atoms.at(j).type));
    if (!li || !lj) throw std::runtime_error("missing LJ parameter in FEP NB14 conversion");
    if (li->first < 0.0 || lj->first < 0.0) throw std::invalid_argument("negative FEP LJ epsilon");
    const double epsilon = std::sqrt(li->first * lj->first);
    const double r6 = lj_combining_rule() == LJCombiningRule::GoodHope
        ? std::pow(4.0 * li->second * lj->second, 3.0)
        : std::pow(li->second + lj->second, 6.0);
    return {epsilon * r6 * r6, 2.0 * epsilon * r6};
}

void exclude_pair(Topology& topology, AtomId i, AtomId j, std::size_t count) {
    if (i >= count || j >= count || i == j) throw std::invalid_argument("invalid FEP explicit pair");
    topology.exclusions.resize(count);
    for (const auto& pair : {std::make_pair(i,j), std::make_pair(j,i)}) {
        auto& row = topology.exclusions[pair.first];
        row.push_back(pair.second);
        std::sort(row.begin(), row.end());
        row.erase(std::unique(row.begin(), row.end()), row.end());
    }
}

std::vector<NB14Extra> endpoint_pairs(const Molecule& molecule) {
    std::vector<NB14Extra> pairs = molecule.nb14_extras;
    for (auto& pair : pairs) {
        if (pair.atom1 >= molecule.atoms.size() || pair.atom2 >= molecule.atoms.size() || pair.atom1 == pair.atom2)
            throw std::invalid_argument("invalid FEP NB14-extra atom index");
        if (!pair.ee_factor) pair.ee_factor = pair.kee * molecule.atoms[pair.atom1].charge * molecule.atoms[pair.atom2].charge;
    }
    for (const auto& pair : build_topology(molecule).nb14s) {
        if (!all_active(molecule, {pair.atom1, pair.atom2})) continue;
        const auto [a,b] = fep_lj_ab(molecule, pair.atom1, pair.atom2);
        pairs.push_back({pair.atom1, pair.atom2, pair.k_lj*a, pair.k_lj*b, pair.k_ee,
                        pair.k_ee*molecule.atoms[pair.atom1].charge*molecule.atoms[pair.atom2].charge});
    }
    return pairs;
}

void merge_explicit_pairs(Molecule& out, const Molecule& a, const Molecule& b, double lambda) {
    struct Coefficients { double a=0., b=0., ee=0., fallback_kee=0.; };
    std::map<std::uint64_t, Coefficients> coefficients;
    const auto append = [&](const Molecule& endpoint, double weight) {
        for (const auto& pair : endpoint_pairs(endpoint)) {
            if (!std::isfinite(pair.a) || !std::isfinite(pair.b) || !std::isfinite(pair.kee) || !std::isfinite(*pair.ee_factor))
                throw std::invalid_argument("non-finite FEP NB14-extra coefficient");
            auto& c = coefficients[pair_key(pair.atom1, pair.atom2)];
            c.a += weight*pair.a; c.b += weight*pair.b;
            c.ee += weight*(*pair.ee_factor); c.fallback_kee += weight*pair.kee;
        }
    };
    append(a,1.0-lambda); append(b,lambda);
    out.nb14_extras.clear();
    out.topology_override->nb14s.clear(); // Converted exactly once, not double counted.
    for (const auto& [key,c] : coefficients) {
        const AtomId i=static_cast<AtomId>(key>>32U), j=static_cast<AtomId>(key&0xffffffffU);
        const double product=out.atoms[i].charge*out.atoms[j].charge;
        // No silent tiny-charge/chemical-state modification. A nonzero physical
        // pair with exactly zero current charge needs another engine representation.
        if (product == 0.0 && c.ee != 0.0)
            throw std::domain_error("FEP NB14 Coulomb coefficient cannot be represented at zero charge product");
        const double kee=product == 0.0 ? c.fallback_kee : c.ee/product;
        if (!std::isfinite(kee)) throw std::domain_error("non-finite FEP NB14 charge compensation");
        out.nb14_extras.push_back({i,j,c.a,c.b,kee,c.ee});
        exclude_pair(*out.topology_override,i,j,out.atoms.size());
    }
}

void require_supported_fep_forces(const Molecule& molecule) {
    if (!molecule.harmonic_impropers.empty() || !molecule.cmaps.empty() || !molecule.urey_bradleys.empty() ||
        !molecule.ryckaert_bellemans.empty() || !molecule.virtual_atoms.empty() || !molecule.listed_force_definitions.empty() ||
        !molecule.sw_parameters.empty() || !molecule.edip_parameters.empty() || molecule.has_gb_parameters)
        throw std::invalid_argument("unsupported force group in FEP merge (ordinary modeling is unaffected)");
}

// Project materialized physical entities, not their union graph. Re-running
// parameter matching after B atoms acquire A/union IDs can choose a different
// ordered improper even when every atom type and parameter is unchanged.
void append_projected_topology(Topology& target, const Topology& source,
                               const std::vector<AtomId>& map,
                               const std::unordered_set<AtomId>* exclusive = nullptr) {
    const auto invalid = std::numeric_limits<AtomId>::max();
    const auto mapped = [&](AtomId id) { return id < map.size() && map[id] != invalid; };
    const auto retain = [&](std::initializer_list<AtomId> ids) {
        return !exclusive || std::any_of(ids.begin(), ids.end(),
            [&](AtomId id) { return exclusive->count(id) != 0; });
    };
    for (auto term : source.bonds) {
        if (!retain({term.atom1, term.atom2})) continue;
        if (!mapped(term.atom1) || !mapped(term.atom2)) continue;
        term.atom1 = map[term.atom1]; term.atom2 = map[term.atom2];
        target.bonds.push_back(term);
    }
    for (auto term : source.angles) {
        if (!retain({term.atom1, term.atom2, term.atom3})) continue;
        if (!mapped(term.atom1) || !mapped(term.atom2) || !mapped(term.atom3)) continue;
        term.atom1 = map[term.atom1]; term.atom2 = map[term.atom2]; term.atom3 = map[term.atom3];
        target.angles.push_back(term);
    }
    for (auto term : source.dihedrals) {
        if (!retain({term.atom1, term.atom2, term.atom3, term.atom4})) continue;
        if (!mapped(term.atom1) || !mapped(term.atom2) || !mapped(term.atom3) || !mapped(term.atom4)) continue;
        term.atom1 = map[term.atom1]; term.atom2 = map[term.atom2];
        term.atom3 = map[term.atom3]; term.atom4 = map[term.atom4];
        target.dihedrals.push_back(term);
    }
    for (auto term : source.nb14s) {
        if (!retain({term.atom1, term.atom2})) continue;
        if (!mapped(term.atom1) || !mapped(term.atom2)) continue;
        term.atom1 = map[term.atom1]; term.atom2 = map[term.atom2];
        target.nb14s.push_back(term);
    }
    for (AtomId id = 0; id < source.exclusions.size(); ++id) {
        if (!mapped(id)) continue;
        for (const auto other : source.exclusions[id])
            if (mapped(other)) target.exclusions[map[id]].push_back(map[other]);
    }
}

void require_disjoint_endpoint(const Topology& topology,
                               const std::function<bool(AtomId)>& inside,
                               const char* message) {
    const auto check = [&](std::initializer_list<AtomId> atoms) {
        const bool first = inside(*atoms.begin());
        for (const auto id : atoms)
            if (inside(id) != first) throw std::invalid_argument(message);
    };
    for (const auto& term : topology.bonds) check({term.atom1, term.atom2});
    for (const auto& term : topology.angles) check({term.atom1, term.atom2, term.atom3});
    for (const auto& term : topology.dihedrals) check({term.atom1, term.atom2, term.atom3, term.atom4});
    for (const auto& term : topology.nb14s) check({term.atom1, term.atom2});
    for (AtomId id = 0; id < topology.exclusions.size(); ++id)
        for (const auto other : topology.exclusions[id]) check({id, other});
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
                    atom.fep_preserve_bonded = true;
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
                atom.fep_preserve_bonded = true;
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
    Topology physical;
    physical.exclusions.resize(out.atoms.size());
    const auto original_a = build_topology(molecule);
    const auto original_b = build_topology(residue_b_molecule);
    std::vector<AtomId> target_map(residue_b_molecule.atoms.size(), invalid);
    std::unordered_set<AtomId> exclusive_a, exclusive_b;
    for (std::uint32_t local = 0; local < residue_a.atom_count; ++local)
        if (a_to_b[local] < 0) exclusive_a.insert(residue_a.atom_begin + local);
    for (std::uint32_t local = 0; local < residue_b.atom_count; ++local) {
        target_map[residue_b.atom_begin + local] = combined_b_to_new[local];
        if (b_to_a[local] < 0) exclusive_b.insert(residue_b.atom_begin + local);
    }
    if (from_a) {
        append_projected_topology(physical, original_a, old_to_new);
        append_projected_topology(physical, original_b, target_map, &exclusive_b);
    } else {
        auto environment_map = old_to_new;
        const auto in_target = [&](AtomId id) {
            return id >= residue_a.atom_begin && id < residue_a.atom_begin + residue_a.atom_count;
        };
        // A standalone replacement residue has no information to reconstruct
        // covalent environment terms. Fail explicitly rather than invent a B
        // Hamiltonian from the union graph; ordinary modeling is unaffected.
        require_disjoint_endpoint(original_a, in_target,
            "covalently linked FEP target requires explicit physical B environment topology");
        for (AtomId id = residue_a.atom_begin; id < residue_a.atom_begin + residue_a.atom_count; ++id)
            environment_map[id] = invalid;
        append_projected_topology(physical, original_a, environment_map);
        require_disjoint_endpoint(original_b, [&](AtomId id) { return target_map[id] != invalid; },
            "covalently linked FEP replacement requires explicit physical B environment topology");
        append_projected_topology(physical, original_b, target_map);
        // Retain exactly the source A branch scaffold, including boundary terms.
        // Do not rebuild from the union graph or weaken any auxiliary coefficient.
        append_projected_topology(physical, original_a, old_to_new, &exclusive_a);
    }
    // Mutually exclusive branches never interact, including through NB14.
    for (const auto a : exclusive_a)
        for (const auto b : exclusive_b)
            exclude_pair(physical, old_to_new[a], target_map[b], out.atoms.size());
    for (auto& excluded : physical.exclusions) {
        std::sort(excluded.begin(), excluded.end());
        excluded.erase(std::unique(excluded.begin(), excluded.end()), excluded.end());
    }
    out.topology_override = std::move(physical);
    return out;
}

}  // namespace

std::pair<Molecule, Molecule> merge_dual_topology(
    const Molecule& molecule, ResidueId residue_index, const Molecule& residue_b_molecule,
    const std::unordered_map<std::uint32_t, std::uint32_t>& match_b_to_a) {
    if (!molecule.validate() || !residue_b_molecule.validate()) {
        throw std::invalid_argument("cannot merge invalid molecule");
    }
    require_supported_fep_forces(molecule);
    require_supported_fep_forces(residue_b_molecule);
    if (!molecule.nb14_extras.empty() || !residue_b_molecule.nb14_extras.empty())
        throw std::invalid_argument("dual topology of pre-existing explicit pairs is not supported; convert after dual assembly");
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
    require_supported_fep_forces(molecule_a);
    require_supported_fep_forces(molecule_b);
    const auto validate_lambda = [](double value) {
        if (!std::isfinite(value) || value < 0.0 || value > 1.0)
            throw std::invalid_argument("FEP lambda must be finite and in [0,1]");
    };
    validate_lambda(default_lambda);
    for (const auto& [name,value] : specific_lambda) {
        if (name!="charge" && name!="mass" && name!="bond_base" && name!="angle" && name!="dihedral")
            throw std::invalid_argument("unsupported FEP specific_lambda key: "+name);
        validate_lambda(value);
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
    out.soft_bonds.clear();
    out.topology_override = merge_topology_for_fep(molecule_a, molecule_b, default_lambda, specific_lambda, out.soft_bonds);
    merge_explicit_pairs(out, molecule_a, molecule_b, lambda_for("dihedral"));

    return out;
}

void fep_intramolecule_nb_to_nb14(Molecule& molecule, const std::vector<ResidueId>& residues) {
    if (!molecule.validate()) throw std::invalid_argument("invalid FEP molecule");
    require_supported_fep_forces(molecule);
    std::set<AtomId> selected;
    for (const auto rid : residues) {
        const auto& residue=molecule.residue(rid); // bounds checked, never silently ignore
        for (AtomId i=residue.atom_begin;i<residue.atom_begin+residue.atom_count;++i) selected.insert(i);
    }
    // Transactional: construct on a copy and replace only after every parameter was found.
    Molecule out=molecule;
    out.topology_override=build_topology(molecule);
    out.topology_override->exclusions.resize(molecule.atoms.size());
    for (auto i=selected.begin();i!=selected.end();++i) {
        for (auto j=std::next(i);j!=selected.end();++j) {
            const auto& x=out.topology_override->exclusions[*i];
            const auto& y=out.topology_override->exclusions[*j];
            const bool explicit_pair=std::any_of(out.nb14_extras.begin(),out.nb14_extras.end(),[&](const NB14Extra& p){return pair_key(p.atom1,p.atom2)==pair_key(*i,*j);});
            if (std::find(x.begin(),x.end(),*j)!=x.end() || std::find(y.begin(),y.end(),*i)!=y.end() || explicit_pair) continue;
            const auto [a,b]=fep_lj_ab(molecule,*i,*j);
            out.nb14_extras.push_back({*i,*j,a,b,1.0,molecule.atoms[*i].charge*molecule.atoms[*j].charge});
            exclude_pair(*out.topology_override,*i,*j,out.atoms.size());
        }
    }
    molecule=std::move(out);
}

Molecule fep_get_free_molecule(const Molecule& molecule, const std::vector<ResidueId>& residues) {
    if (!molecule.validate()) throw std::invalid_argument("invalid FEP molecule");
    require_supported_fep_forces(molecule);
    Molecule out=molecule;
    out.topology_override=build_topology(molecule);
    out.nb14_extras=endpoint_pairs(molecule);
    out.topology_override->nb14s.clear();
    for (const auto rid : residues) {
        const auto& residue=out.residue(rid);
        for (AtomId i=residue.atom_begin;i<residue.atom_begin+residue.atom_count;++i) {
            out.atoms[i].fep_preserve_bonded = !is_zero_lj(out.atoms[i]);
            out.atoms[i].charge=0.0; out.atoms[i].type=kZeroLJ; out.atoms[i].subsys=1;
        }
    }
    // Direct export cannot silently lose stored intramolecular Coulomb energy.
    for (const auto& pair : out.nb14_extras) {
        const double product=out.atoms[pair.atom1].charge*out.atoms[pair.atom2].charge;
        if (product==0.0 && *pair.ee_factor!=0.0)
            throw std::domain_error("decoupled NB14 Coulomb pair at zero charge product is not representable");
    }
    return out;
}

}  // namespace xpongecpp
