#include "core.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>

namespace xpongecpp {
namespace {

using PenaltyTable = std::vector<std::pair<int, int>>;
using BondMap = std::vector<std::unordered_map<std::uint32_t, int>>;
using UnknownConnectivity = std::vector<std::set<std::uint32_t>>;

struct ValencePoint {
    std::uint32_t atom;
    int penalty;
    int valence;
};

struct FormalChargeCheck {
    bool success = false;
    std::vector<std::uint32_t> c3_atoms;
    int choose_c3_negative = -1;
};

const std::map<std::string, PenaltyTable>& atomic_valence() {
    static const std::map<std::string, PenaltyTable> table{
        {"X", {{1, 0}}},
        {"Cn1", {{3, 0}, {4, 1}}},
        {"Cx1", {{4, 0}, {3, 1}}},
        {"C", {{4, 0}, {3, 32}, {2, 64}}},
        {"Nnn1", {{3, 0}, {2, 0}}},
        {"Nx1", {{3, 0}, {2, 3}, {4, 32}}},
        {"Nnn2", {{4, 0}, {3, 1}}},
        {"Nx2", {{3, 0}, {2, 4}, {4, 32}}},
        {"No2", {{5, 0}, {4, 32}, {3, 64}}},
        {"No1", {{5, 0}, {4, 1}, {3, 64}}},
        {"Nx3", {{3, 0}, {4, 1}, {5, 2}, {2, 32}}},
        {"Nx4", {{4, 0}, {3, 64}, {5, 64}}},
        {"On1", {{2, 0}}},
        {"Ox1", {{2, 0}, {1, 1}, {3, 64}}},
        {"Ox2", {{2, 0}, {1, 32}, {3, 64}}},
        {"Px1", {{3, 0}, {2, 2}, {4, 32}}},
        {"Px2", {{3, 0}, {4, 2}, {2, 4}}},
        {"Px3", {{4, 0}, {5, 1}, {6, 2}, {3, 32}}},
        {"Po2", {{5, 0}}},
        {"Po3", {{5, 0}}},
        {"Px4", {{5, 0}, {4, 1}, {6, 32}, {3, 64}}},
        {"Sn1", {{1, 0}, {2, 1}}},
        {"Sx1", {{2, 0}, {1, 2}, {3, 64}}},
        {"Sx2", {{2, 0}, {1, 32}, {3, 32}}},
        {"Sx3", {{4, 0}, {3, 1}, {5, 2}, {6, 2}}},
        {"So2", {{6, 0}}},
        {"So3", {{6, 0}}},
        {"Sx4", {{6, 0}, {5, 2}, {4, 4}}},
        {"Si", {{4, 0}}},
        {"B", {{3, 0}}},
    };
    return table;
}

const std::map<std::string, std::map<int, int>>& atomic_formal_valence() {
    static const std::map<std::string, std::map<int, int>> table{
        {"H", {{1, 0}}},
        {"Cl", {{1, 0}}},
        {"Br", {{1, 0}}},
        {"F", {{1, 0}}},
        {"I", {{1, 0}}},
        {"C", {{4, 0}, {3, 1}, {5, -1}}},
        {"N", {{3, 0}, {4, 1}, {2, -1}, {5, 0}}},
        {"O", {{1, -1}, {2, 0}, {3, 1}}},
        {"P", {{3, 0}, {4, 1}, {2, -1}, {5, 0}, {7, 0}}},
        {"S", {{1, -1}, {2, 0}, {3, 1}, {4, 0}, {6, 0}}},
        {"Si", {{4, 0}}},
        {"B", {{3, 0}}},
    };
    return table;
}

bool bo_atom_judge(const Assign& assign, std::uint32_t atom, const std::string& mask) {
    return assign.atom_judge(atom, mask);
}

int bonded_os_count(const Assign& assign, std::uint32_t atom) {
    int count = 0;
    for (const auto& [neighbor, order] : assign.bonds[atom]) {
        (void)order;
        if (bo_atom_judge(assign, neighbor, "O1") || bo_atom_judge(assign, neighbor, "S1")) {
            ++count;
        }
    }
    return count;
}

std::string determine_bo_type(const Assign& assign, std::uint32_t atom) {
    const auto& element = assign.elements[atom];
    if (element == "H" || element == "F" || element == "Cl" || element == "Br" || element == "I") {
        return "X";
    }
    if (bo_atom_judge(assign, atom, "C1") && assign.elements[assign.bonds[atom].begin()->first] == "N") {
        return "Cn1";
    }
    if (bo_atom_judge(assign, atom, "C1")) {
        return "Cx1";
    }
    if (element == "C") {
        return "C";
    }
    if (bo_atom_judge(assign, atom, "N1") && assign.elements[assign.bonds[atom].begin()->first] == "N") {
        return "Nnn1";
    }
    if (bo_atom_judge(assign, atom, "N1")) {
        return "Nx1";
    }
    if (bo_atom_judge(assign, atom, "N2")) {
        for (const auto& [neighbor, order] : assign.bonds[atom]) {
            (void)order;
            if (assign.elements[neighbor] == "N" && assign.bonds[neighbor].size() == 1) {
                return "Nnn2";
            }
        }
        return "Nx2";
    }
    if (bo_atom_judge(assign, atom, "N3")) {
        const auto count = bonded_os_count(assign, atom);
        if (count == 2) {
            return "No2";
        }
        if (count == 1) {
            return "No1";
        }
        return "Nx3";
    }
    if (bo_atom_judge(assign, atom, "N4")) {
        return "Nx4";
    }
    if (bo_atom_judge(assign, atom, "O1")) {
        int n_count = 0;
        for (const auto& [neighbor, order] : assign.bonds[atom]) {
            (void)order;
            if (assign.elements[neighbor] == "N") {
                ++n_count;
            }
        }
        return n_count == 1 ? "On1" : "Ox1";
    }
    if (bo_atom_judge(assign, atom, "O2")) {
        return "Ox2";
    }
    if (bo_atom_judge(assign, atom, "P1")) {
        return "Px1";
    }
    if (bo_atom_judge(assign, atom, "P2")) {
        return "Px2";
    }
    if (bo_atom_judge(assign, atom, "P3")) {
        return "Px3";
    }
    if (bo_atom_judge(assign, atom, "P4")) {
        const auto count = bonded_os_count(assign, atom);
        if (count == 3) {
            return "Po3";
        }
        if (count == 2) {
            return "Po2";
        }
        return "Px4";
    }
    if (bo_atom_judge(assign, atom, "S1")) {
        int n_count = 0;
        for (const auto& [neighbor, order] : assign.bonds[atom]) {
            (void)order;
            if (assign.elements[neighbor] == "N") {
                ++n_count;
            }
        }
        return n_count == 1 ? "Sn1" : "Sx1";
    }
    if (bo_atom_judge(assign, atom, "S2")) {
        return "Sx2";
    }
    if (bo_atom_judge(assign, atom, "S3")) {
        return "Sx3";
    }
    if (bo_atom_judge(assign, atom, "S4")) {
        const auto count = bonded_os_count(assign, atom);
        if (count >= 3) {
            return "So3";
        }
        if (count == 2) {
            return "So2";
        }
        return "Sx4";
    }
    if (element == "Si") {
        return "Si";
    }
    if (element == "B") {
        return "B";
    }
    throw std::runtime_error("No bond-order atom type for atom #" + std::to_string(atom));
}

std::vector<PenaltyTable> original_penalties(const Assign& assign) {
    std::vector<PenaltyTable> penalties;
    penalties.reserve(assign.elements.size());
    for (std::uint32_t atom = 0; atom < assign.elements.size(); ++atom) {
        const auto type = determine_bo_type(assign, atom);
        penalties.push_back(atomic_valence().at(type));
    }
    return penalties;
}

void collect_unknowns(const Assign& assign, UnknownConnectivity& uc, std::vector<int>& connected) {
    uc.assign(assign.elements.size(), {});
    connected.assign(assign.elements.size(), 0);
    for (std::uint32_t atom = 0; atom < assign.elements.size(); ++atom) {
        for (const auto& [neighbor, order] : assign.bonds[atom]) {
            if (order == -1) {
                uc[atom].insert(neighbor);
            } else if (order > 0) {
                connected[atom] += order;
            }
        }
    }
}

FormalChargeCheck check_formal_charge(Assign& assign, const BondMap& bonds, std::optional<int> total_charge) {
    int total = 0;
    std::vector<std::uint32_t> c3_atoms;
    for (std::uint32_t atom = 0; atom < assign.elements.size(); ++atom) {
        int valence = 0;
        for (const auto& [neighbor, order] : bonds[atom]) {
            (void)neighbor;
            valence += order;
        }
        if (assign.elements[atom] == "C" && valence == 3) {
            c3_atoms.push_back(atom);
        }
        const auto element_it = atomic_formal_valence().find(assign.elements[atom]);
        if (element_it == atomic_formal_valence().end()) {
            return {};
        }
        const auto valence_it = element_it->second.find(valence);
        if (valence_it == element_it->second.end()) {
            return {};
        }
        assign.formal_charges[atom] = valence_it->second;
        total += valence_it->second;
    }
    if (!total_charge || *total_charge == total) {
        return {true, std::move(c3_atoms), -1};
    }
    const int delta_charge = total - *total_charge;
    if (delta_charge % 2 != 0) {
        return {};
    }
    const int need_to_change = delta_charge / 2;
    if (need_to_change == static_cast<int>(c3_atoms.size())) {
        for (const auto atom : c3_atoms) {
            assign.formal_charges[atom] = -1;
        }
        return {true, std::move(c3_atoms), -1};
    }
    if (need_to_change > 0 && need_to_change < static_cast<int>(c3_atoms.size())) {
        return {false, std::move(c3_atoms), need_to_change};
    }
    return {};
}

bool assign_bond_order_one_try(const Assign& assign, const std::vector<int>& valence_input, BondMap& out_bonds) {
    BondMap bonds = assign.bonds;
    UnknownConnectivity uc;
    std::vector<int> connected;
    collect_unknowns(assign, uc, connected);
    std::vector<int> valence = valence_input;
    std::vector<std::vector<int>> valence_backup;
    std::vector<BondMap> bonds_backup;
    std::vector<UnknownConnectivity> uc_backup;
    std::set<std::uint32_t> atom_guessed;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> guess_bonds;
    bool success = false;
    bool determined = false;
    while (!determined && !success) {
        std::vector<std::uint32_t> index_sort(assign.elements.size());
        std::iota(index_sort.begin(), index_sort.end(), 0);
        std::stable_sort(index_sort.begin(), index_sort.end(), [&](std::uint32_t left, std::uint32_t right) {
            const auto left_size = uc[left].empty() ? std::numeric_limits<std::size_t>::max() : uc[left].size();
            const auto right_size = uc[right].empty() ? std::numeric_limits<std::size_t>::max() : uc[right].size();
            return left_size < right_size;
        });
        bool no_basic_rule = true;
        for (const auto i : index_sort) {
            if (uc[i].size() == static_cast<std::size_t>(valence[i]) && valence[i] > 0) {
                while (!uc[i].empty()) {
                    const auto j = *uc[i].rbegin();
                    uc[i].erase(j);
                    bonds[i][j] = 1;
                    bonds[j][i] = 1;
                    uc[j].erase(i);
                    --valence[j];
                }
                valence[i] = 0;
                no_basic_rule = false;
            } else if (uc[i].size() == 1 && valence[i] > 0) {
                const auto j = *uc[i].rbegin();
                uc[i].erase(j);
                uc[j].erase(i);
                bonds[j][i] = valence[i];
                bonds[i][j] = valence[i];
                valence[j] -= valence[i];
                valence[i] = 0;
                no_basic_rule = false;
            }
        }
        success = true;
        for (std::uint32_t i = 0; i < assign.elements.size(); ++i) {
            if (valence[i] != 0 || !uc[i].empty()) {
                success = false;
            }
            if ((uc[i].empty() && valence[i] != 0) || (valence[i] == 0 && !uc[i].empty()) || valence[i] < 0) {
                success = false;
                determined = true;
            }
        }
        if (!success && determined && !guess_bonds.empty()) {
            determined = false;
            const auto [i, j] = guess_bonds.back();
            int trial = bonds[i][j];
            if (trial == 3) {
                guess_bonds.pop_back();
                if (guess_bonds.empty()) {
                    break;
                }
                uc = uc_backup.back();
                uc_backup.pop_back();
                uc[i].insert(j);
                uc[j].insert(i);
                bonds = bonds_backup.back();
                bonds_backup.pop_back();
                valence = valence_backup.back();
                valence_backup.pop_back();
                continue;
            }
            uc = uc_backup.back();
            bonds = bonds_backup.back();
            valence = valence_backup.back();
            ++trial;
            bonds[i][j] = trial;
            bonds[j][i] = trial;
            valence[i] -= trial;
            valence[j] -= trial;
        }
        if (!success && !determined && no_basic_rule) {
            std::vector<std::uint32_t> guess_sort(assign.elements.size());
            std::iota(guess_sort.begin(), guess_sort.end(), 0);
            std::stable_sort(guess_sort.begin(), guess_sort.end(), [&](std::uint32_t left, std::uint32_t right) {
                const int left_value = atom_guessed.count(left) ? -1 : static_cast<int>(uc[left].size());
                const int right_value = atom_guessed.count(right) ? -1 : static_cast<int>(uc[right].size());
                return left_value < right_value;
            });
            const auto i = guess_sort.back();
            if (atom_guessed.count(i) || uc[i].empty()) {
                break;
            }
            const auto j = *uc[i].rbegin();
            uc[i].erase(j);
            atom_guessed.insert(i);
            atom_guessed.insert(j);
            uc[j].erase(i);
            guess_bonds.emplace_back(i, j);
            uc_backup.push_back(uc);
            valence_backup.push_back(valence);
            --valence[i];
            --valence[j];
            bonds_backup.push_back(bonds);
            bonds[i][j] = 1;
            bonds[j][i] = 1;
        }
    }
    if (!success) {
        return false;
    }
    out_bonds = std::move(bonds);
    return true;
}

std::map<int, std::vector<ValencePoint>> collect_penalty_points(const std::vector<PenaltyTable>& penalties) {
    std::map<int, std::vector<ValencePoint>> points;
    for (std::uint32_t atom = 0; atom < penalties.size(); ++atom) {
        for (const auto& [valence, penalty] : penalties[atom]) {
            points[penalty].push_back({atom, penalty, valence});
        }
    }
    return points;
}

std::vector<std::vector<ValencePoint>> preprocess_penalties(
    int penalty,
    const std::map<int, std::vector<ValencePoint>>& penalty_points,
    std::map<int, std::vector<std::vector<ValencePoint>>>& cache) {
    const auto cached = cache.find(penalty);
    if (cached != cache.end()) {
        return cached->second;
    }
    std::vector<std::vector<ValencePoint>> result;
    if (penalty == 1) {
        if (const auto it = penalty_points.find(penalty); it != penalty_points.end()) {
            for (const auto& point : it->second) {
                result.push_back({point});
            }
        }
    } else {
        std::set<std::string> have_added;
        for (int left_penalty = 1; left_penalty <= penalty / 2; ++left_penalty) {
            const auto left = preprocess_penalties(left_penalty, penalty_points, cache);
            const auto right = preprocess_penalties(penalty - left_penalty, penalty_points, cache);
            for (const auto& left_points : left) {
                for (const auto& right_points : right) {
                    std::set<std::uint32_t> atoms;
                    bool duplicate_atom = false;
                    std::vector<ValencePoint> merged = left_points;
                    for (const auto& point : merged) {
                        atoms.insert(point.atom);
                    }
                    for (const auto& point : right_points) {
                        if (atoms.count(point.atom)) {
                            duplicate_atom = true;
                            break;
                        }
                        merged.push_back(point);
                        atoms.insert(point.atom);
                    }
                    if (duplicate_atom) {
                        continue;
                    }
                    std::sort(merged.begin(), merged.end(), [](const ValencePoint& left_point,
                                                               const ValencePoint& right_point) {
                        if (left_point.atom != right_point.atom) {
                            return left_point.atom < right_point.atom;
                        }
                        if (left_point.penalty != right_point.penalty) {
                            return left_point.penalty < right_point.penalty;
                        }
                        return left_point.valence < right_point.valence;
                    });
                    std::string key;
                    for (const auto& point : merged) {
                        if (!key.empty()) {
                            key += "+";
                        }
                        key += std::to_string(point.atom) + "-" + std::to_string(point.penalty) + "-" +
                               std::to_string(point.valence);
                    }
                    if (have_added.insert(key).second) {
                        result.push_back(std::move(merged));
                    }
                }
            }
        }
        if (const auto it = penalty_points.find(penalty); it != penalty_points.end()) {
            for (const auto& point : it->second) {
                result.push_back({point});
            }
        }
    }
    cache[penalty] = result;
    return result;
}

std::vector<std::vector<std::uint32_t>> combinations(const std::vector<std::uint32_t>& values, int choose) {
    std::vector<std::vector<std::uint32_t>> result;
    if (choose < 0 || choose > static_cast<int>(values.size())) {
        return result;
    }
    std::vector<std::uint32_t> current;
    std::function<void(std::size_t, int)> rec = [&](std::size_t start, int remaining) {
        if (remaining == 0) {
            result.push_back(current);
            return;
        }
        for (std::size_t index = start; index <= values.size() - static_cast<std::size_t>(remaining); ++index) {
            current.push_back(values[index]);
            rec(index + 1, remaining - 1);
            current.pop_back();
        }
    };
    rec(0, choose);
    return result;
}

void assign_c3_formal_charge(Assign& assign, const std::vector<std::uint32_t>& c3_atoms,
                             const std::vector<std::uint32_t>& negative_atoms) {
    std::set<std::uint32_t> negative(negative_atoms.begin(), negative_atoms.end());
    for (const auto atom : c3_atoms) {
        assign.formal_charges[atom] = negative.count(atom) ? -1 : 1;
    }
}

}  // namespace

bool Assign::determine_bond_order(bool check_formal_charge_flag, std::optional<int> total_charge) {
    return determine_bond_order_custom(check_formal_charge_flag, total_charge, 2000, 20000, {}, {});
}

bool Assign::determine_bond_order_custom(
    bool check_formal_charge_flag,
    std::optional<int> total_charge,
    int max_step,
    int max_stat,
    const std::vector<std::vector<std::pair<int, int>>>& penalty_scores,
    const std::function<bool(const Assign&)>& extra_criteria) {
    if (max_step <= 0 || max_stat <= 0) {
        throw std::invalid_argument("max_step and max_stat should be positive");
    }
    if (!built) {
        determine_ring_and_bond_type();
    }
    const auto initial_bonds = bonds;
    const auto initial_formal_charges = formal_charges;
    const int input_formal_charge_sum =
        std::accumulate(formal_charges.begin(), formal_charges.end(), 0);
    UnknownConnectivity uc;
    std::vector<int> connected;
    collect_unknowns(*this, uc, connected);
    const auto penalties = penalty_scores.empty() ? original_penalties(*this) : penalty_scores;
    if (penalties.size() != elements.size()) {
        throw std::invalid_argument("penalty_scores should have one entry per atom");
    }
    std::vector<int> best_valence(elements.size(), 0);
    for (std::uint32_t atom = 0; atom < elements.size(); ++atom) {
        if (penalties[atom].empty()) {
            throw std::invalid_argument("penalty_scores entries should not be empty");
        }
        best_valence[atom] = uc[atom].empty() ? 0 : penalties[atom].front().first - connected[atom];
    }
    const auto penalty_points = collect_penalty_points(penalties);
    std::map<int, std::vector<std::vector<ValencePoint>>> cache;

    std::vector<int> valence = best_valence;
    int current_stat = 0;
    std::size_t stat_position = 0;
    std::vector<std::vector<ValencePoint>> points;
    int count = 0;

    BondMap candidate_bonds;
    std::vector<std::vector<std::uint32_t>> formal_charge_combinations;
    std::size_t formal_charge_position = 0;
    std::vector<std::uint32_t> c3_atoms;

    auto get_next_valence = [&]() {
        while (current_stat < max_stat) {
            if (stat_position < points.size()) {
                valence = best_valence;
                bool has_negative_value = false;
                for (const auto& point : points[stat_position]) {
                    valence[point.atom] = point.valence - connected[point.atom];
                    if (valence[point.atom] < 0) {
                        has_negative_value = true;
                        break;
                    }
                }
                ++stat_position;
                if (!has_negative_value) {
                    return true;
                }
            } else {
                ++current_stat;
                stat_position = 0;
                points = preprocess_penalties(current_stat, penalty_points, cache);
            }
        }
        return false;
    };

    auto accept_candidate = [&]() {
        if (!extra_criteria) {
            bonds = std::move(candidate_bonds);
            determine_ring_and_bond_type();
            return true;
        }
        const auto previous_bonds = bonds;
        const auto previous_formal_charges = formal_charges;
        bonds = candidate_bonds;
        if (extra_criteria(*this)) {
            determine_ring_and_bond_type();
            return true;
        }
        bonds = previous_bonds;
        formal_charges = previous_formal_charges;
        return false;
    };

    while (count < max_step && current_stat < max_stat) {
        bool success = false;
        if (!formal_charge_combinations.empty() && formal_charge_position < formal_charge_combinations.size()) {
            assign_c3_formal_charge(*this, c3_atoms, formal_charge_combinations[formal_charge_position]);
            ++formal_charge_position;
            success = true;
        } else {
            formal_charge_combinations.clear();
            formal_charge_position = 0;
            c3_atoms.clear();
            success = assign_bond_order_one_try(*this, valence, candidate_bonds);
            if (success && check_formal_charge_flag) {
                auto formal_check = check_formal_charge(*this, candidate_bonds, total_charge);
                success = formal_check.success;
                if (!success && formal_check.choose_c3_negative > 0) {
                    c3_atoms = std::move(formal_check.c3_atoms);
                    formal_charge_combinations = combinations(c3_atoms, formal_check.choose_c3_negative);
                    if (!formal_charge_combinations.empty()) {
                        continue;
                    }
                }
                if (!success && total_charge.has_value() && !check_connectivity() &&
                    input_formal_charge_sum == *total_charge) {
                    success = true;
                }
            }
        }
        if (success) {
            if (accept_candidate()) {
                return true;
            }
            if (!formal_charge_combinations.empty() &&
                formal_charge_position < formal_charge_combinations.size()) {
                ++count;
                continue;
            }
        }
        ++count;
        if (!get_next_valence()) {
            break;
        }
    }
    bonds = initial_bonds;
    formal_charges = initial_formal_charges;
    return false;
}

}  // namespace xpongecpp
