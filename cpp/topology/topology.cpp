#include "core.hpp"
#include "parallel_utils.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

namespace xpongecpp {
namespace {

struct AtomQuadKey {
    AtomId atom0{0};
    AtomId atom1{0};
    AtomId atom2{0};
    AtomId atom3{0};

    bool operator==(const AtomQuadKey& other) const noexcept {
        return atom0 == other.atom0 && atom1 == other.atom1 &&
               atom2 == other.atom2 && atom3 == other.atom3;
    }
};

struct AtomQuadKeyHash {
    std::size_t operator()(const AtomQuadKey& key) const noexcept {
        std::size_t seed = key.atom0;
        seed = seed * 1315423911U + key.atom1;
        seed = seed * 1315423911U + key.atom2;
        seed = seed * 1315423911U + key.atom3;
        return seed;
    }
};

struct DihedralOccurrence {
    AtomQuadKey key;
    AtomId atom1{0};
    AtomId atom2{0};
    AtomId atom3{0};
    AtomId atom4{0};
};

struct DihedralChunk {
    std::vector<DihedralOccurrence> dihedrals;
    std::vector<std::uint64_t> nb14_pairs;
};

struct PairTypeKey {
    std::uint32_t type1{0};
    std::uint32_t type2{0};

    bool operator==(const PairTypeKey& other) const noexcept {
        return type1 == other.type1 && type2 == other.type2;
    }
};

struct PairTypeKeyHash {
    std::size_t operator()(const PairTypeKey& key) const noexcept {
        return (static_cast<std::size_t>(key.type1) << 32U) ^ key.type2;
    }
};

struct TripleTypeKey {
    std::uint32_t type1{0};
    std::uint32_t type2{0};
    std::uint32_t type3{0};

    bool operator==(const TripleTypeKey& other) const noexcept {
        return type1 == other.type1 && type2 == other.type2 && type3 == other.type3;
    }
};

struct TripleTypeKeyHash {
    std::size_t operator()(const TripleTypeKey& key) const noexcept {
        std::size_t seed = key.type1;
        seed = seed * 1315423911U + key.type2;
        seed = seed * 1315423911U + key.type3;
        return seed;
    }
};

struct QuadTypeKey {
    std::uint32_t type1{0};
    std::uint32_t type2{0};
    std::uint32_t type3{0};
    std::uint32_t type4{0};

    bool operator==(const QuadTypeKey& other) const noexcept {
        return type1 == other.type1 && type2 == other.type2 &&
               type3 == other.type3 && type4 == other.type4;
    }
};

struct QuadTypeKeyHash {
    std::size_t operator()(const QuadTypeKey& key) const noexcept {
        std::size_t seed = key.type1;
        seed = seed * 1315423911U + key.type2;
        seed = seed * 1315423911U + key.type3;
        seed = seed * 1315423911U + key.type4;
        return seed;
    }
};

struct BfsScratch {
    std::vector<std::uint32_t> visit_generation;
    std::vector<std::uint8_t> depth;
    std::vector<AtomId> queue;
    std::uint32_t generation{0};

    void ensure_size(std::size_t atom_count) {
        if (visit_generation.size() != atom_count) {
            visit_generation.assign(atom_count, 0);
            depth.assign(atom_count, 0);
            queue.reserve(atom_count);
            generation = 0;
        }
    }

    void begin(std::size_t atom_count) {
        ensure_size(atom_count);
        ++generation;
        if (generation == 0) {
            std::fill(visit_generation.begin(), visit_generation.end(), 0);
            generation = 1;
        }
        queue.clear();
    }

    bool visited(AtomId atom) const noexcept {
        return visit_generation[atom] == generation;
    }

    void mark(AtomId atom, std::uint8_t atom_depth) {
        visit_generation[atom] = generation;
        depth[atom] = atom_depth;
    }
};

PairTypeKey canonical_pair_key(std::uint32_t type1, std::uint32_t type2) {
    if (type2 < type1) {
        std::swap(type1, type2);
    }
    return {type1, type2};
}

TripleTypeKey canonical_angle_key(std::uint32_t type1, std::uint32_t type2, std::uint32_t type3) {
    if (std::tie(type3, type2, type1) < std::tie(type1, type2, type3)) {
        return {type3, type2, type1};
    }
    return {type1, type2, type3};
}

QuadTypeKey canonical_proper_key(std::uint32_t type1, std::uint32_t type2, std::uint32_t type3, std::uint32_t type4) {
    if (std::tie(type4, type3, type2, type1) < std::tie(type1, type2, type3, type4)) {
        return {type4, type3, type2, type1};
    }
    return {type1, type2, type3, type4};
}

QuadTypeKey ordered_improper_key(std::uint32_t type1, std::uint32_t type2, std::uint32_t type3, std::uint32_t type4) {
    return {type1, type2, type3, type4};
}

AtomQuadKey canonical_improper_group_key(const std::array<AtomId, 4>& atoms) {
    std::array<AtomId, 3> outer{atoms[0], atoms[1], atoms[3]};
    std::sort(outer.begin(), outer.end());
    return {outer[0], outer[1], atoms[2], outer[2]};
}

class TopologyLookupCache {
public:
    explicit TopologyLookupCache(const Molecule& molecule) : molecule_(molecule) {
        std::unordered_map<std::string, std::uint32_t> type_to_id;
        type_to_id.reserve(molecule.atoms.size());
        atom_type_ids_.reserve(molecule.atoms.size());
        std::uint32_t next_id = 1;
        for (const auto& atom : molecule.atoms) {
            const auto [it, inserted] = type_to_id.emplace(atom.type, next_id);
            if (inserted) {
                ++next_id;
            }
            atom_type_ids_.push_back(it->second);
        }
    }

    std::optional<BondTerm> bond(AtomId atom1, AtomId atom2) {
        const auto key = canonical_pair_key(atom_type_ids_[atom1], atom_type_ids_[atom2]);
        {
            std::lock_guard<std::mutex> lock(bond_mutex_);
            const auto it = bond_terms_.find(key);
            if (it != bond_terms_.end()) {
                return it->second;
            }
        }
        const auto result = find_amber_bond_term(molecule_.atoms[atom1].type, molecule_.atoms[atom2].type);
        std::lock_guard<std::mutex> lock(bond_mutex_);
        return bond_terms_.emplace(key, result).first->second;
    }

    std::optional<AngleTerm> angle(AtomId atom1, AtomId atom2, AtomId atom3) {
        const auto key = canonical_angle_key(atom_type_ids_[atom1], atom_type_ids_[atom2], atom_type_ids_[atom3]);
        {
            std::lock_guard<std::mutex> lock(angle_mutex_);
            const auto it = angle_terms_.find(key);
            if (it != angle_terms_.end()) {
                return it->second;
            }
        }
        const auto result = find_amber_angle_term(
            {molecule_.atoms[atom1].type, molecule_.atoms[atom2].type, molecule_.atoms[atom3].type});
        std::lock_guard<std::mutex> lock(angle_mutex_);
        return angle_terms_.emplace(key, result).first->second;
    }

    std::vector<DihedralTerm> proper(AtomId atom1, AtomId atom2, AtomId atom3, AtomId atom4) {
        const auto key = canonical_proper_key(atom_type_ids_[atom1], atom_type_ids_[atom2], atom_type_ids_[atom3],
                                              atom_type_ids_[atom4]);
        {
            std::lock_guard<std::mutex> lock(proper_mutex_);
            const auto it = proper_terms_.find(key);
            if (it != proper_terms_.end()) {
                return it->second;
            }
        }
        const auto result = find_amber_proper_terms(
            {molecule_.atoms[atom1].type, molecule_.atoms[atom2].type, molecule_.atoms[atom3].type, molecule_.atoms[atom4].type});
        std::lock_guard<std::mutex> lock(proper_mutex_);
        return proper_terms_.emplace(key, result).first->second;
    }

    std::optional<AmberImproperMatch> improper(AtomId atom1, AtomId atom2, AtomId atom3, AtomId atom4) {
        const auto key = ordered_improper_key(atom_type_ids_[atom1], atom_type_ids_[atom2], atom_type_ids_[atom3],
                                              atom_type_ids_[atom4]);
        {
            std::lock_guard<std::mutex> lock(improper_mutex_);
            const auto it = improper_terms_.find(key);
            if (it != improper_terms_.end()) {
                return it->second;
            }
        }
        const auto result = find_amber_improper_match(
            {molecule_.atoms[atom1].type, molecule_.atoms[atom2].type, molecule_.atoms[atom3].type, molecule_.atoms[atom4].type});
        std::lock_guard<std::mutex> lock(improper_mutex_);
        return improper_terms_.emplace(key, result).first->second;
    }

    std::optional<NB14Scale> nb14(AtomId atom1, AtomId atom2) {
        const auto key = canonical_pair_key(atom_type_ids_[atom1], atom_type_ids_[atom2]);
        {
            std::lock_guard<std::mutex> lock(nb14_mutex_);
            const auto it = nb14_scales_.find(key);
            if (it != nb14_scales_.end()) {
                return it->second;
            }
        }
        const auto result = find_amber_nb14_scale(molecule_.atoms[atom1].type, molecule_.atoms[atom2].type);
        std::lock_guard<std::mutex> lock(nb14_mutex_);
        return nb14_scales_.emplace(key, result).first->second;
    }

private:
    const Molecule& molecule_;
    std::vector<std::uint32_t> atom_type_ids_;
    std::mutex bond_mutex_;
    std::mutex angle_mutex_;
    std::mutex proper_mutex_;
    std::mutex improper_mutex_;
    std::mutex nb14_mutex_;
    std::unordered_map<PairTypeKey, std::optional<BondTerm>, PairTypeKeyHash> bond_terms_;
    std::unordered_map<TripleTypeKey, std::optional<AngleTerm>, TripleTypeKeyHash> angle_terms_;
    std::unordered_map<QuadTypeKey, std::vector<DihedralTerm>, QuadTypeKeyHash> proper_terms_;
    std::unordered_map<QuadTypeKey, std::optional<AmberImproperMatch>, QuadTypeKeyHash> improper_terms_;
    std::unordered_map<PairTypeKey, std::optional<NB14Scale>, PairTypeKeyHash> nb14_scales_;
};

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
              const Molecule& molecule, TopologyLookupCache* lookup_cache = nullptr) {
    if (atom1 == atom2 || atom1 >= molecule.atoms.size() || atom2 >= molecule.atoms.size()) {
        return;
    }
    const auto key = pair_key(atom1, atom2);
    if (!seen.insert(key).second) {
        return;
    }
    const auto lo = std::min(atom1, atom2);
    const auto hi = std::max(atom1, atom2);
    const auto term = lookup_cache ? lookup_cache->bond(lo, hi)
                                   : find_amber_bond_term(molecule.atoms[lo].type, molecule.atoms[hi].type);
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
    parallel_for_chunks(adjacency.size(), 64, [&](std::size_t begin, std::size_t end) {
        BfsScratch scratch;
        scratch.ensure_size(adjacency.size());
        for (AtomId start = static_cast<AtomId>(begin); start < end; ++start) {
            scratch.begin(adjacency.size());
            scratch.mark(start, 0);
            scratch.queue.push_back(start);
            for (std::size_t head = 0; head < scratch.queue.size(); ++head) {
                const AtomId current = scratch.queue[head];
                if (scratch.depth[current] >= 3) {
                    continue;
                }
                for (const AtomId next : adjacency[current]) {
                    if (scratch.visited(next)) {
                        continue;
                    }
                    const auto next_depth = static_cast<std::uint8_t>(scratch.depth[current] + 1U);
                    scratch.mark(next, next_depth);
                    const int xponge_link_type = next_depth + 1;
                    if (xponge_link_type <= 4) {
                        linked[xponge_link_type][start].push_back(next);
                    }
                    scratch.queue.push_back(next);
                }
            }
            for (auto& layer : linked) {
                std::sort(layer[start].begin(), layer[start].end());
            }
        }
    });
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

AtomQuadKey force_key(const std::array<AtomId, 4>& atoms) {
    return {atoms[0], atoms[1], atoms[2], atoms[3]};
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
                          std::vector<Dihedral>& dihedrals, TopologyLookupCache& lookup_cache) {
    const auto linked = build_linked_atoms(adjacency);
    const std::array<int, 4> topology_like{1, 3, 2, 3};
    const std::array<std::array<int, 4>, 4> topology_matrix{{
        {{1, 3, 2, 3}},
        {{1, 1, 2, 3}},
        {{1, 1, 1, 2}},
        {{1, 1, 1, 1}},
    }};

    const auto candidate_chunks = parallel_collect_chunks(molecule.atoms.size(), 32, [&](std::size_t begin, std::size_t end) {
        std::vector<std::array<AtomId, 4>> local;
        for (AtomId atom0 = static_cast<AtomId>(begin); atom0 < end; ++atom0) {
            std::vector<std::vector<std::array<AtomId, 4>>> backups(4);
            const auto adjacency_degree = adjacency[atom0].size();
            if (adjacency_degree > 0) {
                local.reserve(local.size() + adjacency_degree * adjacency_degree);
            }
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
            local.insert(local.end(), backups[3].begin(), backups[3].end());
        }
        return local;
    });
    std::vector<std::pair<AtomQuadKey, std::array<AtomId, 4>>> grouped_candidates;
    std::size_t total_candidates = 0;
    for (auto& chunk : candidate_chunks) {
        total_candidates += chunk.size();
    }
    grouped_candidates.reserve(total_candidates);
    for (auto& chunk : candidate_chunks) {
        for (const auto& candidate : chunk) {
            grouped_candidates.emplace_back(canonical_improper_group_key(candidate), candidate);
        }
    }
    std::sort(grouped_candidates.begin(), grouped_candidates.end(),
              [](const auto& lhs, const auto& rhs) {
        return std::tie(lhs.first.atom0, lhs.first.atom1, lhs.first.atom2, lhs.first.atom3,
                        lhs.second[0], lhs.second[1], lhs.second[2], lhs.second[3]) <
               std::tie(rhs.first.atom0, rhs.first.atom1, rhs.first.atom2, rhs.first.atom3,
                        rhs.second[0], rhs.second[1], rhs.second[2], rhs.second[3]);
    });

    std::size_t begin = 0;
    while (begin < grouped_candidates.size()) {
        std::size_t end = begin + 1;
        while (end < grouped_candidates.size() && grouped_candidates[end].first == grouped_candidates[begin].first) {
            ++end;
        }
        std::optional<std::pair<std::array<AtomId, 4>, AmberImproperMatch>> selected;
        for (std::size_t index = begin; index < end; ++index) {
            const auto& candidate = grouped_candidates[index].second;
            const auto match = lookup_cache.improper(candidate[0], candidate[1], candidate[2], candidate[3]);
            if (match && match->exact) {
                selected = std::make_pair(candidate, *match);
                break;
            }
        }
        if (!selected) {
            int least_x = 999;
            for (std::size_t index = begin; index < end; ++index) {
                const auto& candidate = grouped_candidates[index].second;
                const auto match = lookup_cache.improper(candidate[0], candidate[1], candidate[2], candidate[3]);
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
        begin = end;
    }
}

}  // namespace

Topology build_topology(const Molecule& molecule) {
    if (!molecule.validate()) {
        throw std::invalid_argument("cannot build topology for invalid molecule");
    }
    if (molecule.topology_override) {
        return *molecule.topology_override;
    }

    Topology topology;
    TopologyLookupCache lookup_cache(molecule);

    std::unordered_set<std::uint64_t> seen_bonds;
    seen_bonds.reserve(molecule.atoms.size() * 2);
    std::vector<bool> residue_has_explicit_bond(molecule.residues.size(), false);
    for (const auto& bond : molecule.explicit_bonds) {
        add_bond(topology.bonds, seen_bonds, bond.atom1, bond.atom2, molecule, &lookup_cache);
        if (bond.atom1 < molecule.atoms.size()) {
            residue_has_explicit_bond[molecule.atoms[bond.atom1].residue] = true;
        }
        if (bond.atom2 < molecule.atoms.size()) {
            residue_has_explicit_bond[molecule.atoms[bond.atom2].residue] = true;
        }
    }
    for (const auto& link : molecule.residue_links) {
        add_bond(topology.bonds, seen_bonds, link.atom1, link.atom2, molecule, &lookup_cache);
    }
    for (ResidueId residue_id = 0; residue_id < molecule.residues.size(); ++residue_id) {
        const auto& residue = molecule.residues[residue_id];
        if (residue_has_explicit_bond[residue_id]) {
            continue;
        }
        if (has_template(residue.name)) {
            const auto& residue_type = get_residue_template(residue.name);
            const auto atoms_by_name = residue_atom_map(molecule, residue);
            if (!molecule.ignore_missing_atoms) {
                for (const auto& template_atom : residue_type.atoms()) {
                    if (atoms_by_name.find(template_atom.name) == atoms_by_name.end()) {
                        throw std::runtime_error(
                            "missing residue atom '" + template_atom.name + "' in residue: " + residue.name);
                    }
                }
            }
            for (const auto& bond : residue_type.bonds()) {
                const auto& atom_name1 = residue_type.atoms()[bond.atom1].name;
                const auto& atom_name2 = residue_type.atoms()[bond.atom2].name;
                const auto it1 = atoms_by_name.find(atom_name1);
                const auto it2 = atoms_by_name.find(atom_name2);
                if (it1 != atoms_by_name.end() && it2 != atoms_by_name.end()) {
                    add_bond(topology.bonds, seen_bonds, it1->second, it2->second, molecule, &lookup_cache);
                }
            }
            continue;
        }
        if (residue.name == "NA" || residue.name == "CL") {
            continue;
        }
        if (residue.atom_count > 1) {
            throw std::runtime_error("missing residue template/connectivity for residue: " + residue.name);
        }
    }

    std::sort(topology.bonds.begin(), topology.bonds.end(), [](const Bond& lhs, const Bond& rhs) {
        return std::tie(lhs.atom1, lhs.atom2) < std::tie(rhs.atom1, rhs.atom2);
    });

    const auto adjacency = build_adjacency(molecule.atoms.size(), topology.bonds);

    const auto angle_chunks = parallel_collect_chunks(adjacency.size(), 32, [&](std::size_t begin, std::size_t end) {
        std::vector<Angle> local;
        for (AtomId center = static_cast<AtomId>(begin); center < end; ++center) {
            if (molecule.residues[molecule.atoms[center].residue].name == "WAT") {
                continue;
            }
            const auto& neighbors = adjacency[center];
            if (neighbors.size() < 2) {
                continue;
            }
            local.reserve(local.size() + (neighbors.size() * (neighbors.size() - 1)) / 2);
            for (std::size_t i = 0; i < neighbors.size(); ++i) {
                for (std::size_t k = i + 1; k < neighbors.size(); ++k) {
                    const auto term = lookup_cache.angle(neighbors[i], center, neighbors[k]);
                    if (term && term->k != 0.0) {
                        local.push_back({neighbors[i], center, neighbors[k], term->k, term->theta});
                    } else if (!term) {
                        local.push_back({neighbors[i], center, neighbors[k], 40.0, 1.91});
                    }
                }
            }
        }
        return local;
    });
    for (auto& chunk : angle_chunks) {
        topology.angles.insert(topology.angles.end(), chunk.begin(), chunk.end());
    }

    const auto dihedral_chunks = parallel_collect_chunks(topology.bonds.size(), 16, [&](std::size_t begin, std::size_t end) {
        DihedralChunk local;
        for (std::size_t bond_index = begin; bond_index < end; ++bond_index) {
            const auto& bond = topology.bonds[bond_index];
            const AtomId j = bond.atom1;
            const AtomId k = bond.atom2;
            local.dihedrals.reserve(local.dihedrals.size() + adjacency[j].size() * adjacency[k].size());
            for (const AtomId i : adjacency[j]) {
                if (i == k) {
                    continue;
                }
                for (const AtomId l : adjacency[k]) {
                    if (l == j || l == i) {
                        continue;
                    }
                    const bool swap_outer = l < i;
                    local.dihedrals.push_back({
                        {swap_outer ? l : i, j, k, swap_outer ? i : l},
                        i,
                        j,
                        k,
                        l,
                    });
                    if (!has_shorter_path_than_three(adjacency, i, l)) {
                        local.nb14_pairs.push_back(pair_key(i, l));
                    }
                }
            }
        }
        std::sort(local.dihedrals.begin(), local.dihedrals.end(),
                  [](const DihedralOccurrence& lhs, const DihedralOccurrence& rhs) {
            return std::tie(lhs.key.atom0, lhs.key.atom1, lhs.key.atom2, lhs.key.atom3,
                            lhs.atom1, lhs.atom2, lhs.atom3, lhs.atom4) <
                   std::tie(rhs.key.atom0, rhs.key.atom1, rhs.key.atom2, rhs.key.atom3,
                            rhs.atom1, rhs.atom2, rhs.atom3, rhs.atom4);
        });
        local.dihedrals.erase(
            std::unique(local.dihedrals.begin(), local.dihedrals.end(),
                        [](const DihedralOccurrence& lhs, const DihedralOccurrence& rhs) {
                return lhs.key == rhs.key &&
                       lhs.atom1 == rhs.atom1 && lhs.atom2 == rhs.atom2 &&
                       lhs.atom3 == rhs.atom3 && lhs.atom4 == rhs.atom4;
            }),
            local.dihedrals.end());
        std::sort(local.nb14_pairs.begin(), local.nb14_pairs.end());
        local.nb14_pairs.erase(std::unique(local.nb14_pairs.begin(), local.nb14_pairs.end()), local.nb14_pairs.end());
        return local;
    });

    std::vector<DihedralOccurrence> merged_dihedrals;
    std::size_t total_dihedral_occurrences = 0;
    std::size_t total_nb14_pairs = 0;
    for (const auto& chunk : dihedral_chunks) {
        total_dihedral_occurrences += chunk.dihedrals.size();
        total_nb14_pairs += chunk.nb14_pairs.size();
    }
    merged_dihedrals.reserve(total_dihedral_occurrences);
    std::vector<std::uint64_t> seen_14;
    seen_14.reserve(total_nb14_pairs);
    for (const auto& chunk : dihedral_chunks) {
        merged_dihedrals.insert(merged_dihedrals.end(), chunk.dihedrals.begin(), chunk.dihedrals.end());
        seen_14.insert(seen_14.end(), chunk.nb14_pairs.begin(), chunk.nb14_pairs.end());
    }
    std::sort(merged_dihedrals.begin(), merged_dihedrals.end(),
              [](const DihedralOccurrence& lhs, const DihedralOccurrence& rhs) {
        return std::tie(lhs.key.atom0, lhs.key.atom1, lhs.key.atom2, lhs.key.atom3,
                        lhs.atom1, lhs.atom2, lhs.atom3, lhs.atom4) <
               std::tie(rhs.key.atom0, rhs.key.atom1, rhs.key.atom2, rhs.key.atom3,
                        rhs.atom1, rhs.atom2, rhs.atom3, rhs.atom4);
    });
    merged_dihedrals.erase(
        std::unique(merged_dihedrals.begin(), merged_dihedrals.end(),
                    [](const DihedralOccurrence& lhs, const DihedralOccurrence& rhs) {
            return lhs.key == rhs.key;
        }),
        merged_dihedrals.end());
    for (const auto& occurrence : merged_dihedrals) {
        const auto terms = lookup_cache.proper(occurrence.atom1, occurrence.atom2, occurrence.atom3, occurrence.atom4);
        for (const auto& term : terms) {
            push_oriented_dihedral(topology.dihedrals, occurrence.atom1, occurrence.atom2,
                                   occurrence.atom3, occurrence.atom4, term);
        }
    }
    std::sort(seen_14.begin(), seen_14.end());
    seen_14.erase(std::unique(seen_14.begin(), seen_14.end()), seen_14.end());

    add_xponge_impropers(molecule, adjacency, topology.dihedrals, lookup_cache);

    std::sort(topology.dihedrals.begin(), topology.dihedrals.end(), [](const Dihedral& lhs, const Dihedral& rhs) {
        return std::tie(lhs.atom1, lhs.atom2, lhs.atom3, lhs.atom4, lhs.periodicity, lhs.k, lhs.phase) <
               std::tie(rhs.atom1, rhs.atom2, rhs.atom3, rhs.atom4, rhs.periodicity, rhs.k, rhs.phase);
    });

    topology.exclusions.resize(molecule.atoms.size());
    parallel_for_chunks(adjacency.size(), 32, [&](std::size_t begin, std::size_t end) {
        BfsScratch scratch;
        scratch.ensure_size(adjacency.size());
        for (AtomId start = static_cast<AtomId>(begin); start < end; ++start) {
            auto& excluded = topology.exclusions[start];
            excluded.clear();
            scratch.begin(adjacency.size());
            scratch.mark(start, 0);
            scratch.queue.push_back(start);
            for (std::size_t head = 0; head < scratch.queue.size(); ++head) {
                const AtomId current = scratch.queue[head];
                if (scratch.depth[current] >= 3) {
                    continue;
                }
                for (const AtomId next : adjacency[current]) {
                    if (scratch.visited(next)) {
                        continue;
                    }
                    const auto next_depth = static_cast<std::uint8_t>(scratch.depth[current] + 1U);
                    scratch.mark(next, next_depth);
                    scratch.queue.push_back(next);
                    excluded.push_back(next);
                }
            }
            std::sort(excluded.begin(), excluded.end());
        }
    });

    topology.nb14s.reserve(seen_14.size());
    for (const auto key : seen_14) {
        const AtomId atom1 = static_cast<AtomId>(key >> 32U);
        const AtomId atom2 = static_cast<AtomId>(key & 0xffffffffU);
        if (const auto scale = lookup_cache.nb14(atom1, atom2)) {
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
