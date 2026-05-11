#include "core.hpp"

#include <algorithm>
#include <cmath>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

namespace xpongecpp {
namespace {

std::uint64_t pair_key(AtomId atom1, AtomId atom2) {
    const auto lo = std::min(atom1, atom2);
    const auto hi = std::max(atom1, atom2);
    return (static_cast<std::uint64_t>(lo) << 32U) | hi;
}

double distance(const Atom& atom1, const Atom& atom2) {
    const double dx = atom1.x - atom2.x;
    const double dy = atom1.y - atom2.y;
    const double dz = atom1.z - atom2.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

bool is_solvent_or_ion(const Residue& residue) {
    return residue.name == "WAT" || residue.name == "HOH" || residue.name == "NA" || residue.name == "CL";
}

bool likely_bonded(const Atom& atom1, const Atom& atom2) {
    if (atom1.element == "H" && atom2.element == "H") {
        return false;
    }
    const double cutoff = (atom1.element == "H" || atom2.element == "H") ? 1.25 : 1.95;
    return distance(atom1, atom2) <= cutoff;
}

AtomId find_atom(const Molecule& molecule, const Residue& residue, const std::string& name) {
    for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
        const AtomId atom_id = residue.atom_begin + local;
        if (molecule.atoms[atom_id].name == name) {
            return atom_id;
        }
    }
    return static_cast<AtomId>(molecule.atoms.size());
}

std::unordered_map<std::string, AtomId> residue_atom_map(const Molecule& molecule, const Residue& residue) {
    std::unordered_map<std::string, AtomId> out;
    out.reserve(residue.atom_count);
    for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
        const AtomId atom_id = residue.atom_begin + local;
        out[molecule.atoms[atom_id].name] = atom_id;
    }
    return out;
}

void add_bond(std::vector<Bond>& bonds, std::unordered_set<std::uint64_t>& seen, AtomId atom1, AtomId atom2,
              const Molecule& molecule) {
    if (atom1 == atom2 || atom1 >= molecule.atoms.size() || atom2 >= molecule.atoms.size()) {
        return;
    }
    const auto key = pair_key(atom1, atom2);
    if (!seen.insert(key).second) {
        return;
    }
    const auto lo = std::min(atom1, atom2);
    const auto hi = std::max(atom1, atom2);
    const auto term = find_amber_bond_term(molecule.atoms[lo].type, molecule.atoms[hi].type);
    if (term) {
        bonds.push_back({lo, hi, term->k, term->length});
    } else {
        bonds.push_back({lo, hi, 300.0, distance(molecule.atoms[lo], molecule.atoms[hi])});
    }
}

std::vector<std::vector<AtomId>> build_adjacency(std::size_t atom_count, const std::vector<Bond>& bonds) {
    std::vector<std::vector<AtomId>> adjacency(atom_count);
    for (const auto& bond : bonds) {
        adjacency[bond.atom1].push_back(bond.atom2);
        adjacency[bond.atom2].push_back(bond.atom1);
    }
    for (auto& neighbors : adjacency) {
        std::sort(neighbors.begin(), neighbors.end());
    }
    return adjacency;
}

bool has_shorter_path_than_three(const std::vector<std::vector<AtomId>>& adjacency, AtomId atom1, AtomId atom2) {
    if (std::binary_search(adjacency[atom1].begin(), adjacency[atom1].end(), atom2)) {
        return true;
    }
    for (const auto neighbor : adjacency[atom1]) {
        if (std::binary_search(adjacency[atom2].begin(), adjacency[atom2].end(), neighbor)) {
            return true;
        }
    }
    return false;
}

std::array<std::vector<std::vector<AtomId>>, 5> build_linked_atoms(const std::vector<std::vector<AtomId>>& adjacency) {
    std::array<std::vector<std::vector<AtomId>>, 5> linked;
    for (auto& layer : linked) {
        layer.resize(adjacency.size());
    }
    for (AtomId start = 0; start < adjacency.size(); ++start) {
        std::vector<int> depth(adjacency.size(), -1);
        std::queue<AtomId> queue;
        depth[start] = 0;
        queue.push(start);
        while (!queue.empty()) {
            const AtomId current = queue.front();
            queue.pop();
            if (depth[current] >= 3) {
                continue;
            }
            for (const AtomId next : adjacency[current]) {
                if (depth[next] != -1) {
                    continue;
                }
                depth[next] = depth[current] + 1;
                const int xponge_link_type = depth[next] + 1;
                if (xponge_link_type <= 4) {
                    linked[xponge_link_type][start].push_back(next);
                }
                queue.push(next);
            }
        }
        for (auto& layer : linked) {
            std::sort(layer[start].begin(), layer[start].end());
        }
    }
    return linked;
}

bool contains_atom(const std::vector<AtomId>& atoms, AtomId atom) {
    return std::binary_search(atoms.begin(), atoms.end(), atom);
}

std::array<std::string, 4> atom_types_for(const Molecule& molecule, AtomId atom1, AtomId atom2, AtomId atom3,
                                          AtomId atom4) {
    return {molecule.atoms[atom1].type, molecule.atoms[atom2].type, molecule.atoms[atom3].type,
            molecule.atoms[atom4].type};
}

void push_oriented_dihedral(std::vector<Dihedral>& dihedrals, AtomId atom1, AtomId atom2, AtomId atom3, AtomId atom4,
                            const DihedralTerm& term) {
    if (term.k == 0.0) {
        return;
    }
    if (atom1 > atom4) {
        std::swap(atom1, atom4);
        std::swap(atom2, atom3);
    }
    dihedrals.push_back({atom1, atom2, atom3, atom4, term.periodicity, term.k, term.phase});
}

std::string force_key(const std::array<AtomId, 4>& atoms) {
    return std::to_string(atoms[0]) + ":" + std::to_string(atoms[1]) + ":" + std::to_string(atoms[2]) + ":" +
           std::to_string(atoms[3]);
}

std::vector<std::array<AtomId, 4>> improper_same_force(const std::array<AtomId, 4>& atoms) {
    std::array<AtomId, 3> outer{atoms[0], atoms[1], atoms[3]};
    std::sort(outer.begin(), outer.end());
    std::vector<std::array<AtomId, 4>> out;
    do {
        out.push_back({outer[0], outer[1], atoms[2], outer[2]});
    } while (std::next_permutation(outer.begin(), outer.end()));
    return out;
}

void add_xponge_impropers(const Molecule& molecule, const std::vector<std::vector<AtomId>>& adjacency,
                          std::vector<Dihedral>& dihedrals) {
    const auto linked = build_linked_atoms(adjacency);
    const std::array<int, 4> topology_like{1, 3, 2, 3};
    const std::array<std::array<int, 4>, 4> topology_matrix{{
        {{1, 3, 2, 3}},
        {{1, 1, 2, 3}},
        {{1, 1, 1, 2}},
        {{1, 1, 1, 1}},
    }};

    std::vector<std::array<AtomId, 4>> candidates;
    for (AtomId atom0 = 0; atom0 < molecule.atoms.size(); ++atom0) {
        std::vector<std::vector<std::array<AtomId, 4>>> backups(4);
        backups[0].push_back({atom0, 0, 0, 0});
        for (std::size_t i = 1; i < topology_like.size(); ++i) {
            for (const AtomId atom1 : linked[topology_like[i]][atom0]) {
                for (const auto& backup : backups[i - 1]) {
                    bool good = true;
                    for (std::size_t j = 0; j < i; ++j) {
                        const auto distance = std::abs(topology_matrix[j][i]);
                        if (backup[j] == atom1 || distance <= 1 || !contains_atom(linked[distance][backup[j]], atom1)) {
                            good = false;
                            break;
                        }
                        if (topology_matrix[j][i] <= -1) {
                            for (int shorter = 2; shorter < distance; ++shorter) {
                                if (contains_atom(linked[shorter][backup[j]], atom1)) {
                                    good = false;
                                    break;
                                }
                            }
                        }
                    }
                    if (good) {
                        auto next = backup;
                        next[i] = atom1;
                        backups[i].push_back(next);
                    }
                }
            }
        }
        candidates.insert(candidates.end(), backups[3].begin(), backups[3].end());
    }

    std::vector<std::vector<std::array<AtomId, 4>>> groups;
    std::unordered_map<std::string, std::size_t> group_by_key;
    for (const auto& candidate : candidates) {
        const auto key = force_key(candidate);
        const auto it = group_by_key.find(key);
        if (it != group_by_key.end()) {
            groups[it->second].push_back(candidate);
            continue;
        }
        const auto group_index = groups.size();
        groups.push_back({candidate});
        for (const auto& equivalent : improper_same_force(candidate)) {
            group_by_key[force_key(equivalent)] = group_index;
        }
    }

    for (const auto& group : groups) {
        std::optional<std::pair<std::array<AtomId, 4>, AmberImproperMatch>> selected;
        for (const auto& candidate : group) {
            const auto match = find_amber_improper_match(atom_types_for(molecule, candidate[0], candidate[1],
                                                                        candidate[2], candidate[3]));
            if (match && match->exact) {
                selected = std::make_pair(candidate, *match);
                break;
            }
        }
        if (!selected) {
            int least_x = 999;
            for (const auto& candidate : group) {
                const auto match = find_amber_improper_match(atom_types_for(molecule, candidate[0], candidate[1],
                                                                            candidate[2], candidate[3]));
                if (!match || match->exact || match->wildcard_count > least_x) {
                    continue;
                }
                selected = std::make_pair(candidate, *match);
                least_x = match->wildcard_count;
            }
        }
        if (selected && selected->second.term.k != 0.0) {
            const auto& atoms = selected->first;
            push_oriented_dihedral(dihedrals, atoms[0], atoms[1], atoms[2], atoms[3], selected->second.term);
        }
    }
}

}  // namespace

Topology build_topology(const Molecule& molecule) {
    Topology topology;
    if (!molecule.validate()) {
        throw std::invalid_argument("cannot build topology for invalid molecule");
    }

    std::unordered_set<std::uint64_t> seen_bonds;
    seen_bonds.reserve(molecule.atoms.size() * 2);
    std::vector<bool> residue_has_explicit_bond(molecule.residues.size(), false);
    for (const auto& bond : molecule.explicit_bonds) {
        add_bond(topology.bonds, seen_bonds, bond.atom1, bond.atom2, molecule);
        if (bond.atom1 < molecule.atoms.size()) {
            residue_has_explicit_bond[molecule.atoms[bond.atom1].residue] = true;
        }
        if (bond.atom2 < molecule.atoms.size()) {
            residue_has_explicit_bond[molecule.atoms[bond.atom2].residue] = true;
        }
    }
    for (ResidueId residue_id = 0; residue_id < molecule.residues.size(); ++residue_id) {
        const auto& residue = molecule.residues[residue_id];
        if (residue_has_explicit_bond[residue_id]) {
            continue;
        }
        if (has_template(residue.name)) {
            const auto& residue_type = get_residue_template(residue.name);
            const auto atoms_by_name = residue_atom_map(molecule, residue);
            for (const auto& bond : residue_type.bonds()) {
                const auto& atom_name1 = residue_type.atoms()[bond.atom1].name;
                const auto& atom_name2 = residue_type.atoms()[bond.atom2].name;
                const auto it1 = atoms_by_name.find(atom_name1);
                const auto it2 = atoms_by_name.find(atom_name2);
                if (it1 != atoms_by_name.end() && it2 != atoms_by_name.end()) {
                    add_bond(topology.bonds, seen_bonds, it1->second, it2->second, molecule);
                }
            }
            continue;
        }
        if (residue.name == "NA" || residue.name == "CL") {
            continue;
        }
        for (std::uint32_t local1 = 0; local1 < residue.atom_count; ++local1) {
            const AtomId atom1 = residue.atom_begin + local1;
            for (std::uint32_t local2 = local1 + 1; local2 < residue.atom_count; ++local2) {
                const AtomId atom2 = residue.atom_begin + local2;
                if (likely_bonded(molecule.atoms[atom1], molecule.atoms[atom2])) {
                    add_bond(topology.bonds, seen_bonds, atom1, atom2, molecule);
                }
            }
        }
    }

    for (std::size_t residue_index = 0; residue_index + 1 < molecule.residues.size(); ++residue_index) {
        const auto& current = molecule.residues[residue_index];
        const auto& next = molecule.residues[residue_index + 1];
        if (residue_has_explicit_bond[residue_index] || residue_has_explicit_bond[residue_index + 1]) {
            continue;
        }
        if (is_solvent_or_ion(current) || is_solvent_or_ion(next)) {
            continue;
        }
        const AtomId carbon = find_atom(molecule, current, "C");
        const AtomId nitrogen = find_atom(molecule, next, "N");
        if (carbon < molecule.atoms.size() && nitrogen < molecule.atoms.size()) {
            add_bond(topology.bonds, seen_bonds, carbon, nitrogen, molecule);
        }
    }

    std::sort(topology.bonds.begin(), topology.bonds.end(), [](const Bond& lhs, const Bond& rhs) {
        return std::tie(lhs.atom1, lhs.atom2) < std::tie(rhs.atom1, rhs.atom2);
    });

    const auto adjacency = build_adjacency(molecule.atoms.size(), topology.bonds);

    for (AtomId center = 0; center < adjacency.size(); ++center) {
        if (molecule.residues[molecule.atoms[center].residue].name == "WAT") {
            continue;
        }
        const auto& neighbors = adjacency[center];
        for (std::size_t i = 0; i < neighbors.size(); ++i) {
            for (std::size_t k = i + 1; k < neighbors.size(); ++k) {
                const auto term = find_amber_angle_term(
                    {molecule.atoms[neighbors[i]].type, molecule.atoms[center].type, molecule.atoms[neighbors[k]].type});
                if (term && term->k != 0.0) {
                    topology.angles.push_back({neighbors[i], center, neighbors[k], term->k, term->theta});
                } else if (!term) {
                    topology.angles.push_back({neighbors[i], center, neighbors[k], 40.0, 1.91});
                }
            }
        }
    }

    std::unordered_set<std::uint64_t> seen_14;
    std::unordered_set<std::string> seen_dihedrals;
    for (const auto& bond : topology.bonds) {
        const AtomId j = bond.atom1;
        const AtomId k = bond.atom2;
        for (const AtomId i : adjacency[j]) {
            if (i == k) {
                continue;
            }
            for (const AtomId l : adjacency[k]) {
                if (l == j || l == i) {
                    continue;
                }
                const auto a = std::min(i, l);
                const auto d = std::max(i, l);
                const std::string key = std::to_string(a) + ":" + std::to_string(j) + ":" +
                                        std::to_string(k) + ":" + std::to_string(d);
                if (seen_dihedrals.insert(key).second) {
                    const auto terms = find_amber_proper_terms(atom_types_for(molecule, i, j, k, l));
                    for (const auto& term : terms) {
                        push_oriented_dihedral(topology.dihedrals, i, j, k, l, term);
                    }
                }
                if (!has_shorter_path_than_three(adjacency, i, l)) {
                    seen_14.insert(pair_key(i, l));
                }
            }
        }
    }

    add_xponge_impropers(molecule, adjacency, topology.dihedrals);

    std::sort(topology.dihedrals.begin(), topology.dihedrals.end(), [](const Dihedral& lhs, const Dihedral& rhs) {
        return std::tie(lhs.atom1, lhs.atom2, lhs.atom3, lhs.atom4, lhs.periodicity, lhs.k, lhs.phase) <
               std::tie(rhs.atom1, rhs.atom2, rhs.atom3, rhs.atom4, rhs.periodicity, rhs.k, rhs.phase);
    });

    topology.exclusions.resize(molecule.atoms.size());
    for (AtomId start = 0; start < adjacency.size(); ++start) {
        std::vector<int> depth(adjacency.size(), -1);
        std::queue<AtomId> queue;
        depth[start] = 0;
        queue.push(start);
        while (!queue.empty()) {
            const AtomId current = queue.front();
            queue.pop();
            if (depth[current] >= 3) {
                continue;
            }
            for (const AtomId next : adjacency[current]) {
                if (depth[next] != -1) {
                    continue;
                }
                depth[next] = depth[current] + 1;
                queue.push(next);
                topology.exclusions[start].push_back(next);
            }
        }
        std::sort(topology.exclusions[start].begin(), topology.exclusions[start].end());
    }

    topology.nb14s.reserve(seen_14.size());
    for (const auto key : seen_14) {
        const AtomId atom1 = static_cast<AtomId>(key >> 32U);
        const AtomId atom2 = static_cast<AtomId>(key & 0xffffffffU);
        if (const auto scale = find_amber_nb14_scale(molecule.atoms[atom1].type, molecule.atoms[atom2].type)) {
            if (scale->k_lj != 0.0 && scale->k_ee != 0.0) {
                topology.nb14s.push_back({atom1, atom2, scale->k_lj, scale->k_ee});
            }
        }
    }
    std::sort(topology.nb14s.begin(), topology.nb14s.end(), [](const NB14& lhs, const NB14& rhs) {
        return std::tie(lhs.atom1, lhs.atom2) < std::tie(rhs.atom1, rhs.atom2);
    });

    return topology;
}

}  // namespace xpongecpp
