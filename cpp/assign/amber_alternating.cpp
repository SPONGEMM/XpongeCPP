#include "amber_alternating.hpp"

#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace xpongecpp {
namespace {

using BondTriple = std::tuple<std::uint32_t, std::uint32_t, int>;

std::vector<BondTriple> iter_sequence_bonds(const Assign& assign) {
    std::vector<BondTriple> bonds;
    std::set<std::pair<std::uint32_t, std::uint32_t>> seen;
    for (const auto& pair : assign.bond_sequence) {
        const auto atom_i = pair.first;
        const auto atom_j = pair.second;
        if (atom_i >= assign.bonds.size()) {
            continue;
        }
        const auto found = assign.bonds[atom_i].find(atom_j);
        if (found == assign.bonds[atom_i].end()) {
            continue;
        }
        seen.insert({atom_i, atom_j});
        bonds.emplace_back(atom_i, atom_j, found->second);
    }
    for (std::uint32_t atom_i = 0; atom_i < assign.bonds.size(); ++atom_i) {
        for (const auto& [atom_j, bond_order] : assign.bonds[atom_i]) {
            if (atom_i < atom_j && seen.count({atom_i, atom_j}) == 0) {
                bonds.emplace_back(atom_i, atom_j, bond_order);
            }
        }
    }
    return bonds;
}

void normalize_to_primary(Assign& assign, const std::unordered_map<std::string, std::string>& primary_to_secondary) {
    std::unordered_map<std::string, std::string> secondary_to_primary;
    for (const auto& [primary, secondary] : primary_to_secondary) {
        secondary_to_primary[secondary] = primary;
    }
    for (auto& atom_type : assign.atom_types) {
        const auto found = secondary_to_primary.find(atom_type);
        if (found != secondary_to_primary.end()) {
            atom_type = found->second;
        }
    }
}

void orient_and_convert(Assign& assign,
                        const std::unordered_map<std::string, std::string>& primary_to_secondary,
                        const std::unordered_set<int>& same_bond_orders,
                        const std::unordered_set<int>& opposite_bond_orders) {
    normalize_to_primary(assign, primary_to_secondary);
    std::unordered_set<std::string> candidate_names;
    for (const auto& [primary, secondary] : primary_to_secondary) {
        candidate_names.insert(primary);
    }
    std::vector<int> atom_sign(assign.atom_count(), 0);
    std::vector<bool> is_candidate(assign.atom_count(), false);
    bool has_candidate = false;
    int num = 0;
    for (std::uint32_t atom_i = 0; atom_i < assign.atom_count(); ++atom_i) {
        if (candidate_names.count(assign.atom_types[atom_i]) == 0) {
            continue;
        }
        is_candidate[atom_i] = true;
        if (!has_candidate) {
            atom_sign[atom_i] = 1;
            has_candidate = true;
        }
        ++num;
    }
    if (!has_candidate) {
        return;
    }
    --num;
    const auto sequence_bonds = iter_sequence_bonds(assign);
    while (num > 0) {
        --num;
        int flag = 0;
        for (const auto& [atom_i, atom_j, bond_order] : sequence_bonds) {
            if (!(is_candidate[atom_i] && is_candidate[atom_j])) {
                continue;
            }
            if (flag == 0 && atom_sign[atom_i] == 0 && atom_sign[atom_j] == 0) {
                atom_sign[atom_i] = 1;
            }
            if (atom_sign[atom_i] == 0 && atom_sign[atom_j] != 0) {
                flag = 1;
                if (same_bond_orders.count(bond_order) != 0) {
                    atom_sign[atom_i] = atom_sign[atom_j];
                } else if (opposite_bond_orders.count(bond_order) != 0) {
                    atom_sign[atom_i] = -atom_sign[atom_j];
                }
            }
            if (atom_sign[atom_j] == 0 && atom_sign[atom_i] != 0) {
                flag = 1;
                if (same_bond_orders.count(bond_order) != 0) {
                    atom_sign[atom_j] = atom_sign[atom_i];
                } else if (opposite_bond_orders.count(bond_order) != 0) {
                    atom_sign[atom_j] = -atom_sign[atom_i];
                }
            }
        }
    }

    for (std::uint32_t atom_i = 0; atom_i < assign.atom_count(); ++atom_i) {
        if (atom_sign[atom_i] != -1) {
            continue;
        }
        const auto found = primary_to_secondary.find(assign.atom_types[atom_i]);
        if (found != primary_to_secondary.end()) {
            assign.atom_types[atom_i] = found->second;
        }
    }
}

void adjust_cp_cq(Assign& assign) {
    normalize_to_primary(assign, {{"cp", "cq"}});
    std::vector<int> atom_sign(assign.atom_count(), 0);
    std::vector<bool> is_candidate(assign.atom_count(), false);
    int index = 0;
    int num = 0;
    for (std::uint32_t atom_i = 0; atom_i < assign.atom_count(); ++atom_i) {
        if (assign.atom_types[atom_i] != "cp") {
            continue;
        }
        is_candidate[atom_i] = true;
        if (index == 0) {
            atom_sign[atom_i] = 1;
            index = 1;
        }
        ++num;
    }
    if (num == 0) {
        return;
    }
    --num;
    const auto sequence_bonds = iter_sequence_bonds(assign);
    while (num > 0) {
        --num;
        for (const auto& [atom_i, atom_j, bond_order] : sequence_bonds) {
            if (!(is_candidate[atom_i] && is_candidate[atom_j])) {
                continue;
            }
            const bool same_phase = bond_order == 1 && !assign.has_bond_marker(atom_i, atom_j, "AB");
            if (atom_sign[atom_i] == 0 && atom_sign[atom_j] != 0) {
                atom_sign[atom_i] = same_phase ? atom_sign[atom_j] : -atom_sign[atom_j];
            }
            if (atom_sign[atom_j] == 0 && atom_sign[atom_i] != 0) {
                atom_sign[atom_j] = same_phase ? atom_sign[atom_i] : -atom_sign[atom_i];
            }
        }
    }
    for (std::uint32_t atom_i = 0; atom_i < assign.atom_count(); ++atom_i) {
        if (atom_sign[atom_i] == -1 && assign.atom_types[atom_i] == "cp") {
            assign.atom_types[atom_i] = "cq";
        }
    }
}

void demote_specific_nc_nd_to_n2(Assign& assign) {
    for (std::uint32_t atom_i = 0; atom_i < assign.atom_count(); ++atom_i) {
        if (assign.atom_types[atom_i] != "nc" && assign.atom_types[atom_i] != "nd") {
            continue;
        }
        if (!assign.has_atom_marker(atom_i, "AR3")) {
            continue;
        }
        if (assign.atom_marker_count(atom_i, "RG5") != 1 || assign.has_atom_marker(atom_i, "RG9")) {
            continue;
        }
        bool has_single_n3_neighbor = false;
        for (const auto& [bonded_atom, bond_order] : assign.bonds[atom_i]) {
            if (bond_order == 1 && assign.atom_judge(bonded_atom, "N3")) {
                has_single_n3_neighbor = true;
                break;
            }
        }
        if (has_single_n3_neighbor) {
            assign.atom_types[atom_i] = "n2";
        }
    }
}

void demote_sequence_sensitive_ne_nf_to_n2(Assign& assign) {
    std::set<std::pair<std::uint32_t, std::uint32_t>> sequence_bonds(
        assign.bond_sequence.begin(), assign.bond_sequence.end()
    );
    for (std::uint32_t atom_i = 0; atom_i < assign.atom_count(); ++atom_i) {
        if (assign.atom_types[atom_i] != "ne" && assign.atom_types[atom_i] != "nf") {
            continue;
        }
        std::vector<std::uint32_t> terminal_double_heteros;
        bool has_other_terminal_hetero = false;
        for (const auto& [bonded_atom, bond_order] : assign.bonds[atom_i]) {
            const auto& element = assign.elements[bonded_atom];
            if ((element != "O" && element != "S") || assign.bonds[bonded_atom].size() != 1) {
                continue;
            }
            if (bond_order == 2) {
                terminal_double_heteros.push_back(bonded_atom);
            } else {
                has_other_terminal_hetero = true;
            }
        }
        if (terminal_double_heteros.size() != 1 || has_other_terminal_hetero) {
            continue;
        }
        const auto hetero = terminal_double_heteros.front();
        if (sequence_bonds.count({hetero, atom_i}) != 0) {
            assign.atom_types[atom_i] = "n2";
        }
    }
}

}  // namespace

void apply_amber_alternating_type_adjustment(Assign& assign) {
    demote_specific_nc_nd_to_n2(assign);
    demote_sequence_sensitive_ne_nf_to_n2(assign);
    orient_and_convert(
        assign,
        {
            {"cc", "cd"},
            {"ce", "cf"},
            {"cg", "ch"},
            {"pc", "pd"},
            {"pe", "pf"},
            {"nc", "nd"},
            {"ne", "nf"},
        },
        {1, 7},
        {2, 3, 8}
    );
    adjust_cp_cq(assign);
}

}  // namespace xpongecpp
