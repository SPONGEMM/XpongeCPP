#include "core.hpp"

#include <algorithm>
#include <numeric>
#include <set>

namespace xpongecpp {
namespace {

struct RingInfo {
    std::vector<std::uint32_t> atoms;
    bool is_pure_aromatic_ring{false};
    bool is_pure_aliphatic_ring{false};
    bool is_planar_ring{false};
    bool out_plane_double_bond{false};
};

std::vector<std::uint32_t> canonical_ring(std::vector<std::uint32_t> ring) {
    const auto min_it = std::min_element(ring.begin(), ring.end());
    std::rotate(ring.begin(), min_it, ring.end());
    auto reversed = ring;
    std::reverse(reversed.begin(), reversed.end());
    std::rotate(reversed.begin(), reversed.end() - 1, reversed.end());
    if (reversed.size() > 1 && reversed[1] < ring[1]) {
        ring = std::move(reversed);
    }
    return ring;
}

std::string ring_key(const std::vector<std::uint32_t>& ring) {
    std::string key;
    for (const auto atom : ring) {
        if (!key.empty()) {
            key.push_back('-');
        }
        key += std::to_string(atom);
    }
    return key;
}

std::vector<RingInfo> get_rings(const Assign& assign) {
    std::set<std::string> seen;
    std::vector<RingInfo> rings;
    for (std::uint32_t atom0 = 0; atom0 < assign.elements.size(); ++atom0) {
        std::vector<std::uint32_t> current_path{atom0};
        std::unordered_map<std::uint32_t, int> current_path_sons;
        std::vector<std::pair<std::uint32_t, std::uint32_t>> current_work;
        for (const auto& [atom, order] : assign.bonds[atom0]) {
            (void)order;
            current_work.emplace_back(atom, atom0);
        }
        current_path_sons[atom0] = static_cast<int>(assign.bonds[atom0].size());
        std::vector<std::uint32_t> current_path_father;
        while (!current_path.empty() && !current_work.empty()) {
            const auto [work_atom, from_atom] = current_work.back();
            current_work.pop_back();
            current_path.push_back(work_atom);
            current_path_father.push_back(from_atom);

            std::vector<std::pair<std::uint32_t, std::uint32_t>> bond_atom;
            for (const auto& [atom, order] : assign.bonds[work_atom]) {
                (void)order;
                if (atom == from_atom) {
                    continue;
                }
                const auto found = std::find(current_path.begin(), current_path.end(), atom);
                if (found != current_path.end()) {
                    std::vector<std::uint32_t> ring(found, current_path.end());
                    ring = canonical_ring(std::move(ring));
                    const auto key = ring_key(ring);
                    if (seen.insert(key).second) {
                        rings.push_back({std::move(ring)});
                    }
                } else {
                    bond_atom.emplace_back(atom, work_atom);
                }
            }

            if (current_path.size() < 9) {
                current_path_sons[work_atom] = static_cast<int>(bond_atom.size());
                current_work.insert(current_work.end(), bond_atom.begin(), bond_atom.end());
            } else {
                current_path_sons[work_atom] = 0;
            }

            for (auto it = current_path.rbegin(); it != current_path.rend();) {
                const auto atom = *it;
                if (current_path_sons[atom] != 0) {
                    break;
                }
                const auto pop_atom = current_path.back();
                current_path.pop_back();
                current_path_sons.erase(pop_atom);
                if (!current_path_father.empty()) {
                    const auto father_atom = current_path_father.back();
                    current_path_father.pop_back();
                    --current_path_sons[father_atom];
                }
                it = current_path.rbegin();
            }
        }
    }
    return rings;
}

void check_pure_aromatic(const Assign& assign, RingInfo& ring) {
    if (ring.atoms.size() != 6) {
        ring.is_pure_aromatic_ring = false;
        return;
    }
    ring.is_pure_aromatic_ring = true;
    for (const auto atom : ring.atoms) {
        if (!assign.atom_judge(atom, "C3") && !assign.atom_judge(atom, "N2") && !assign.atom_judge(atom, "N3")) {
            ring.is_pure_aromatic_ring = false;
            break;
        }
        if (assign.atom_judge(atom, "N3")) {
            int total_order = 0;
            for (const auto& [bonded_atom, bond_order] : assign.bonds[atom]) {
                (void)bonded_atom;
                total_order += bond_order;
            }
            if (total_order == 3) {
                ring.is_pure_aromatic_ring = false;
                break;
            }
        }
        for (const auto& [bonded_atom, bond_order] : assign.bonds[atom]) {
            if (bond_order == 2 && !assign.has_atom_marker(bonded_atom, "RG")) {
                ring.is_pure_aromatic_ring = false;
                break;
            }
        }
        if (!ring.is_pure_aromatic_ring) {
            break;
        }
    }
}

void check_pure_aliphatic_and_planar(Assign& assign, RingInfo& ring) {
    ring.is_pure_aliphatic_ring = true;
    ring.is_planar_ring = true;
    for (const auto atom : ring.atoms) {
        if (ring.is_pure_aromatic_ring) {
            assign.add_atom_marker(atom, "AR1");
            for (std::size_t i = 0; i < 6; ++i) {
                const auto atom1 = ring.atoms[(i + 5) % 6];
                const auto atom2 = ring.atoms[i];
                assign.add_bond_marker(atom1, atom2, "AB");
            }
        }
        if (!assign.atom_judge(atom, "C4")) {
            ring.is_pure_aliphatic_ring = false;
        }
        if (!assign.atom_judge(atom, "C3") && !assign.atom_judge(atom, "N2") && !assign.atom_judge(atom, "N3") &&
            !assign.atom_judge(atom, "O2") && !assign.atom_judge(atom, "S2") && !assign.atom_judge(atom, "P2") &&
            !assign.atom_judge(atom, "P3")) {
            ring.is_planar_ring = false;
        }
    }
}

void check_out_plane_double_bond(const Assign& assign, RingInfo& ring) {
    ring.out_plane_double_bond = false;
    for (const auto atom : ring.atoms) {
        for (const auto& [bonded_atom, order] : assign.bonds[atom]) {
            if (assign.elements[bonded_atom] != "C" && order == 2 &&
                std::find(ring.atoms.begin(), ring.atoms.end(), bonded_atom) == ring.atoms.end()) {
                ring.out_plane_double_bond = true;
            }
        }
    }
}

}  // namespace

void Assign::determine_ring_and_bond_type() {
    (void)check_connectivity();
    for (std::uint32_t atom = 0; atom < elements.size(); ++atom) {
        atom_markers[atom].clear();
        int dlo = 0;
        int noto = 0;
        for (const auto& [atom2, order] : bonds[atom]) {
            (void)order;
            bond_markers[atom][atom2].clear();
            if (atom_judge(atom2, "O1")) {
                ++dlo;
            } else {
                ++noto;
            }
        }
        if (dlo >= 1 && noto <= 1) {
            for (const auto& [atom2, order] : bonds[atom]) {
                (void)order;
                if (atom_judge(atom2, "O1")) {
                    add_bond_marker(atom, atom2, "DLB");
                }
            }
        }
        for (const auto& [atom2, order] : bonds[atom]) {
            if (has_bond_marker(atom, atom2, "DLB")) {
                add_bond_marker(atom, atom2, "DL", true);
                add_bond_marker(atom, atom2, "sb", true);
            } else if (order == 1) {
                add_bond_marker(atom, atom2, "sb", true);
                if (!has_bond_marker(atom, atom2, "AB")) {
                    add_bond_marker(atom, atom2, "SB", true);
                }
            } else if (order == 2) {
                add_bond_marker(atom, atom2, "db", true);
                if (!has_bond_marker(atom, atom2, "AB")) {
                    add_bond_marker(atom, atom2, "DB", true);
                }
            } else {
                add_bond_marker(atom, atom2, "tb", true);
            }
        }
    }

    auto ring_infos = get_rings(*this);
    rings.clear();
    rings.reserve(ring_infos.size());
    for (const auto& ring : ring_infos) {
        rings.push_back(ring.atoms);
        for (const auto atom : ring.atoms) {
            add_atom_marker(atom, "RG");
            add_atom_marker(atom, "RG" + std::to_string(ring.atoms.size()));
        }
    }
    for (auto& ring : ring_infos) {
        check_pure_aromatic(*this, ring);
        check_pure_aliphatic_and_planar(*this, ring);
        check_out_plane_double_bond(*this, ring);
        if (!ring.is_pure_aromatic_ring) {
            for (const auto atom : ring.atoms) {
                if (ring.is_pure_aliphatic_ring) {
                    add_atom_marker(atom, "AR5");
                } else if (ring.is_planar_ring) {
                    add_atom_marker(atom, ring.out_plane_double_bond ? "AR3" : "AR2");
                } else {
                    add_atom_marker(atom, "AR4");
                }
            }
        }
    }
    built = true;
}

}  // namespace xpongecpp
