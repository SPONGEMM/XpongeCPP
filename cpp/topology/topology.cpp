#include "core.hpp"

#include <algorithm>
#include <cmath>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
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
    bonds.push_back({lo, hi, 300.0, distance(molecule.atoms[lo], molecule.atoms[hi])});
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

}  // namespace

Topology build_topology(const Molecule& molecule) {
    Topology topology;
    if (!molecule.validate()) {
        throw std::invalid_argument("cannot build topology for invalid molecule");
    }

    std::unordered_set<std::uint64_t> seen_bonds;
    seen_bonds.reserve(molecule.atoms.size() * 2);
    for (const auto& residue : molecule.residues) {
        if (residue.name == "WAT" || residue.name == "HOH") {
            const AtomId oxygen = find_atom(molecule, residue, "O");
            const AtomId h1 = find_atom(molecule, residue, "H1");
            const AtomId h2 = find_atom(molecule, residue, "H2");
            if (oxygen < molecule.atoms.size() && h1 < molecule.atoms.size()) {
                add_bond(topology.bonds, seen_bonds, oxygen, h1, molecule);
            }
            if (oxygen < molecule.atoms.size() && h2 < molecule.atoms.size()) {
                add_bond(topology.bonds, seen_bonds, oxygen, h2, molecule);
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
        const auto& neighbors = adjacency[center];
        for (std::size_t i = 0; i < neighbors.size(); ++i) {
            for (std::size_t k = i + 1; k < neighbors.size(); ++k) {
                topology.angles.push_back({neighbors[i], center, neighbors[k], 40.0, 1.91});
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
                    topology.dihedrals.push_back({i, j, k, l, 3, 0.1, 3.141592653589793});
                }
                seen_14.insert(pair_key(i, l));
            }
        }
    }

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
        topology.nb14s.push_back({atom1, atom2, 0.5, 0.833333});
    }
    std::sort(topology.nb14s.begin(), topology.nb14s.end(), [](const NB14& lhs, const NB14& rhs) {
        return std::tie(lhs.atom1, lhs.atom2) < std::tie(rhs.atom1, rhs.atom2);
    });

    return topology;
}

}  // namespace xpongecpp
