#include "core.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace xpongecpp {
namespace {

constexpr double kAngstromPerBohr = 0.52918;
constexpr double kPi = 3.14159265358979323846;

double norm3(const std::array<double, 3>& value) {
    return std::sqrt(value[0] * value[0] + value[1] * value[1] + value[2] * value[2]);
}

std::array<double, 3> subtract3(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

std::vector<std::array<double, 3>> fibonacci_grid(int npoints, const std::array<double, 3>& center, double radius) {
    std::vector<std::array<double, 3>> out;
    if (npoints <= 0) {
        return out;
    }
    out.reserve(static_cast<std::size_t>(npoints));
    const double golden_angle = kPi * (3.0 - std::sqrt(5.0));
    for (int i = 0; i < npoints; ++i) {
        const double y = 1.0 - (2.0 * i + 1.0) / static_cast<double>(npoints);
        const double r = std::sqrt(std::max(1.0 - y * y, 0.0));
        const double theta = golden_angle * static_cast<double>(i);
        out.push_back({
            center[0] + radius * std::cos(theta) * r,
            center[1] + radius * y,
            center[2] + radius * std::sin(theta) * r,
        });
    }
    return out;
}

std::vector<double> solve_linear_system(std::vector<std::vector<double>> matrix, std::vector<double> rhs) {
    const std::size_t n = matrix.size();
    if (n == 0 || rhs.size() != n) {
        throw std::invalid_argument("linear system dimensions do not match");
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (matrix[i].size() != n) {
            throw std::invalid_argument("linear system matrix must be square");
        }
        matrix[i].push_back(rhs[i]);
    }
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        double pivot_value = std::abs(matrix[pivot][col]);
        for (std::size_t row = col + 1; row < n; ++row) {
            const double candidate = std::abs(matrix[row][col]);
            if (candidate > pivot_value) {
                pivot = row;
                pivot_value = candidate;
            }
        }
        if (pivot_value < 1e-14) {
            throw std::runtime_error("singular RESP linear system");
        }
        if (pivot != col) {
            std::swap(matrix[pivot], matrix[col]);
        }
        const double diag = matrix[col][col];
        for (std::size_t j = col; j <= n; ++j) {
            matrix[col][j] /= diag;
        }
        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) {
                continue;
            }
            const double factor = matrix[row][col];
            if (factor == 0.0) {
                continue;
            }
            for (std::size_t j = col; j <= n; ++j) {
                matrix[row][j] -= factor * matrix[col][j];
            }
        }
    }
    std::vector<double> solution(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        solution[i] = matrix[i][n];
    }
    return solution;
}

std::vector<double> force_equivalence_q(std::vector<double> q, const std::vector<std::vector<int>>& extra_equivalence) {
    for (const auto& group : extra_equivalence) {
        if (group.empty()) {
            continue;
        }
        double mean = 0.0;
        for (int atom : group) {
            mean += q.at(static_cast<std::size_t>(atom));
        }
        mean /= static_cast<double>(group.size());
        for (int atom : group) {
            q.at(static_cast<std::size_t>(atom)) = mean;
        }
    }
    return q;
}

bool atom_judge(const Assign& assign, std::uint32_t atom, const std::string& mask) {
    std::string element;
    std::string digits;
    for (char ch : mask) {
        if (std::isdigit(static_cast<unsigned char>(ch))) {
            digits.push_back(ch);
        } else {
            element.push_back(ch);
        }
    }
    if (!digits.empty()) {
        return assign.elements.at(atom) == element &&
               assign.bonds.at(atom).size() == static_cast<std::size_t>(std::stoi(digits));
    }
    return assign.elements.at(atom) == element;
}

void find_tofit_second(
    const Assign& assign,
    std::vector<std::vector<int>>& tofit_second,
    std::vector<int>& fit_group,
    int& sublength
) {
    const int atom_count = static_cast<int>(assign.atom_count());
    fit_group.assign(static_cast<std::size_t>(atom_count), -1);
    sublength = 0;
    for (int i = 0; i < atom_count; ++i) {
        if (atom_judge(assign, static_cast<std::uint32_t>(i), "C4")) {
            fit_group[static_cast<std::size_t>(i)] = static_cast<int>(tofit_second.size());
            tofit_second.push_back({i});
            std::vector<int> hydrogens;
            for (const auto& bonded : assign.bonds[static_cast<std::size_t>(i)]) {
                if (assign.elements[bonded.first] == "H") {
                    hydrogens.push_back(static_cast<int>(bonded.first));
                }
            }
            if (!hydrogens.empty()) {
                for (int atom : hydrogens) {
                    fit_group[static_cast<std::size_t>(atom)] = static_cast<int>(tofit_second.size());
                }
                tofit_second.push_back(hydrogens);
                sublength += static_cast<int>(hydrogens.size()) - 1;
            }
        }
        if (atom_judge(assign, static_cast<std::uint32_t>(i), "C3")) {
            std::vector<int> hydrogens;
            for (const auto& bonded : assign.bonds[static_cast<std::size_t>(i)]) {
                if (assign.elements[bonded.first] == "H") {
                    hydrogens.push_back(static_cast<int>(bonded.first));
                }
            }
            if (hydrogens.size() == 2) {
                fit_group[static_cast<std::size_t>(i)] = static_cast<int>(tofit_second.size());
                tofit_second.push_back({i});
                for (int atom : hydrogens) {
                    fit_group[static_cast<std::size_t>(atom)] = static_cast<int>(tofit_second.size());
                }
                tofit_second.push_back(hydrogens);
                sublength += 1;
            }
        }
    }
}

void correct_extra_equivalence(
    std::vector<std::vector<int>>& tofit_second,
    std::vector<int>& fit_group,
    int& sublength,
    const std::vector<std::vector<int>>& extra_equivalence
) {
    if (extra_equivalence.empty()) {
        return;
    }
    std::vector<std::vector<int>> equi_group;
    equi_group.reserve(extra_equivalence.size());
    for (const auto& eq : extra_equivalence) {
        std::vector<int> group;
        for (int atom : eq) {
            const int value = fit_group.at(static_cast<std::size_t>(atom));
            if (value != -1 && std::find(group.begin(), group.end(), value) == group.end()) {
                group.push_back(value);
            }
        }
        std::sort(group.begin(), group.end());
        equi_group.push_back(group);
    }
    std::vector<int> all_groups = fit_group;
    std::sort(all_groups.begin(), all_groups.end());
    all_groups.erase(std::unique(all_groups.begin(), all_groups.end()), all_groups.end());
    std::unordered_map<int, int> group_map;
    for (int group : all_groups) {
        group_map[group] = group;
    }
    for (const auto& eq : equi_group) {
        if (eq.empty()) {
            continue;
        }
        for (int group : eq) {
            group_map[group] = eq.front();
        }
    }
    int temp_max = 0;
    for (int group : all_groups) {
        if (group == -1) {
            continue;
        }
        if (group_map[group] == group) {
            group_map[group] = temp_max++;
        } else {
            group_map[group] = group_map[group_map[group]];
        }
    }
    auto old = tofit_second;
    tofit_second.assign(static_cast<std::size_t>(temp_max), {});
    for (std::size_t i = 0; i < old.size(); ++i) {
        auto& target = tofit_second[static_cast<std::size_t>(group_map[static_cast<int>(i)])];
        target.insert(target.end(), old[i].begin(), old[i].end());
        sublength -= static_cast<int>(old[i].size()) - 1;
    }
    for (const auto& group : tofit_second) {
        sublength += static_cast<int>(group.size()) - 1;
    }
    for (std::size_t atom = 0; atom < fit_group.size(); ++atom) {
        fit_group[atom] = group_map[fit_group[atom]];
    }
}

std::vector<int> restrained_second_stage_groups(
    const Assign& assign,
    const std::vector<std::vector<int>>& tofit_second
) {
    std::vector<int> restrained;
    for (std::size_t group_index = 0; group_index < tofit_second.size(); ++group_index) {
        const auto& group = tofit_second[group_index];
        const bool has_heavy_atom = std::any_of(group.begin(), group.end(), [&assign](int atom) {
            return assign.elements.at(static_cast<std::size_t>(atom)) != "H";
        });
        if (has_heavy_atom) {
            restrained.push_back(static_cast<int>(group_index));
        }
    }
    return restrained;
}

std::vector<double> resp_scf_kernel(
    const Assign& assign,
    int atom_count,
    double a,
    double b,
    std::vector<std::vector<double>> matrix_a,
    const std::vector<std::vector<double>>& matrix_a0,
    const std::vector<double>& matrix_b,
    std::vector<double> q
) {
    std::vector<double> q_last_step = q;
    int step = 0;
    double max_delta = 0.0;
    do {
        ++step;
        q_last_step = q;
        for (int i = 0; i < atom_count; ++i) {
            if (assign.elements[static_cast<std::size_t>(i)] != "H") {
                matrix_a[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] =
                    matrix_a0[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] +
                    a / std::sqrt(q_last_step[static_cast<std::size_t>(i)] * q_last_step[static_cast<std::size_t>(i)] + b * b);
            }
        }
        auto solution = solve_linear_system(matrix_a, matrix_b);
        q.assign(solution.begin(), solution.end() - 1);
        max_delta = 0.0;
        for (std::size_t i = 0; i < q.size(); ++i) {
            max_delta = std::max(max_delta, std::abs(q[i] - q_last_step[i]));
        }
    } while (step == 1 || max_delta > 1e-4);
    return q;
}

void get_a20_and_b20(
    int total_length,
    const std::vector<std::vector<int>>& tofit_second,
    std::vector<int> fit_group,
    int sublength,
    int atom_count,
    const std::vector<std::vector<double>>& matrix_a0,
    const std::vector<double>& matrix_b,
    int charge,
    const std::vector<double>& q,
    std::vector<std::vector<double>>& a20,
    std::vector<double>& b20
) {
    a20.assign(static_cast<std::size_t>(total_length), std::vector<double>(static_cast<std::size_t>(total_length), 0.0));
    int count = static_cast<int>(tofit_second.size());
    for (int i = 0; i < atom_count; ++i) {
        if (fit_group[static_cast<std::size_t>(i)] == -1) {
            fit_group[static_cast<std::size_t>(i)] = count++;
        }
        a20[static_cast<std::size_t>(atom_count - sublength)][static_cast<std::size_t>(fit_group[static_cast<std::size_t>(i)])] += 1.0;
        a20[static_cast<std::size_t>(fit_group[static_cast<std::size_t>(i)])][static_cast<std::size_t>(atom_count - sublength)] += 1.0;
    }
    b20.assign(static_cast<std::size_t>(total_length), 0.0);
    for (int i = 0; i < atom_count; ++i) {
        b20[static_cast<std::size_t>(fit_group[static_cast<std::size_t>(i)])] += matrix_b[static_cast<std::size_t>(i)];
        for (int j = 0; j < atom_count; ++j) {
            a20[static_cast<std::size_t>(fit_group[static_cast<std::size_t>(i)])][static_cast<std::size_t>(fit_group[static_cast<std::size_t>(j)])] +=
                matrix_a0[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        }
    }
    b20[static_cast<std::size_t>(atom_count - sublength)] = static_cast<double>(charge);
    count = 0;
    for (int i = 0; i < atom_count; ++i) {
        if (fit_group[static_cast<std::size_t>(i)] >= static_cast<int>(tofit_second.size())) {
            b20[static_cast<std::size_t>(atom_count - sublength + count + 1)] = q[static_cast<std::size_t>(i)];
            a20[static_cast<std::size_t>(atom_count - sublength + count + 1)][static_cast<std::size_t>(tofit_second.size() + count)] = 1.0;
            a20[static_cast<std::size_t>(tofit_second.size() + count)][static_cast<std::size_t>(atom_count - sublength + count + 1)] = 1.0;
            ++count;
        }
    }
}

}  // namespace

std::vector<std::array<double, 3>> generate_resp_mk_grid(
    const std::vector<std::string>& atoms,
    const std::vector<std::array<double, 3>>& atom_coordinates_bohr,
    double area_density,
    int layer,
    const std::unordered_map<std::string, double>& radius
) {
    static const std::unordered_map<std::string, double> default_radius = {
        {"H", 1.2}, {"C", 1.5}, {"N", 1.5}, {"O", 1.4}, {"P", 1.8}, {"S", 1.75},
        {"F", 1.35}, {"Cl", 1.7}, {"Br", 2.3},
    };
    std::unordered_map<std::string, double> real_radius = default_radius;
    for (const auto& item : radius) {
        real_radius[item.first] = item.second;
    }
    std::vector<std::array<double, 3>> grids;
    const double factor = area_density * kAngstromPerBohr * kAngstromPerBohr * 4.0 * kPi;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        const auto found = real_radius.find(atoms[i]);
        if (found == real_radius.end()) {
            throw std::out_of_range("Radius for element " + atoms[i] + " not found");
        }
        const double r0 = found->second / kAngstromPerBohr;
        for (int shell = 0; shell < layer; ++shell) {
            const double r = r0 * (1.4 + 0.2 * static_cast<double>(shell));
            auto shell_points = fibonacci_grid(static_cast<int>(factor * r * r), atom_coordinates_bohr[i], r);
            grids.insert(grids.end(), shell_points.begin(), shell_points.end());
        }
    }
    std::vector<std::array<double, 3>> filtered;
    filtered.reserve(grids.size());
    for (const auto& point : grids) {
        bool keep = true;
        for (std::size_t i = 0; i < atoms.size(); ++i) {
            const double r0 = 1.39 * real_radius.at(atoms[i]) / kAngstromPerBohr;
            if (norm3(subtract3(point, atom_coordinates_bohr[i])) < r0) {
                keep = false;
                break;
            }
        }
        if (keep) {
            filtered.push_back(point);
        }
    }
    return filtered;
}

RespFitDebugResult fit_resp_from_esp_cpp_debug(
    const Assign& assign,
    const std::vector<std::array<double, 3>>& atom_coordinates_bohr,
    const std::vector<double>& nuclear_charges,
    const std::vector<std::array<double, 3>>& grid_points_bohr,
    const std::vector<double>& esp_values_au,
    int charge,
    const std::vector<std::vector<int>>& extra_equivalence,
    double a1,
    double a2,
    bool two_stage,
    bool only_esp
) {
    const auto total_start = std::chrono::steady_clock::now();
    const int atom_count = static_cast<int>(assign.atom_count());
    RespFitDebugResult result;
    const auto assembly_start = std::chrono::steady_clock::now();
    std::vector<std::vector<double>> matrix_a0(static_cast<std::size_t>(atom_count), std::vector<double>(static_cast<std::size_t>(atom_count), 0.0));
    std::vector<double> vnuc(grid_points_bohr.size(), 0.0);
    for (int i = 0; i < atom_count; ++i) {
        for (std::size_t grid = 0; grid < grid_points_bohr.size(); ++grid) {
            const auto rp = subtract3(atom_coordinates_bohr[static_cast<std::size_t>(i)], grid_points_bohr[grid]);
            const double norm_rp = norm3(rp);
            vnuc[grid] += nuclear_charges[static_cast<std::size_t>(i)] / norm_rp;
        }
        for (int j = 0; j < atom_count; ++j) {
            double value = 0.0;
            for (std::size_t grid = 0; grid < grid_points_bohr.size(); ++grid) {
                const double norm_i = norm3(subtract3(atom_coordinates_bohr[static_cast<std::size_t>(i)], grid_points_bohr[grid]));
                const double norm_j = norm3(subtract3(atom_coordinates_bohr[static_cast<std::size_t>(j)], grid_points_bohr[grid]));
                value += 1.0 / norm_i / norm_j;
            }
            matrix_a0[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = value;
        }
    }
    for (auto& row : matrix_a0) {
        row.push_back(1.0);
    }
    matrix_a0.push_back(std::vector<double>(static_cast<std::size_t>(atom_count + 1), 1.0));
    matrix_a0.back().back() = 0.0;
    auto matrix_a = matrix_a0;

    std::vector<double> mep(grid_points_bohr.size(), 0.0);
    for (std::size_t i = 0; i < grid_points_bohr.size(); ++i) {
        mep[i] = vnuc[i] - esp_values_au[i];
    }
    std::vector<double> matrix_b(static_cast<std::size_t>(atom_count + 1), 0.0);
    for (int i = 0; i < atom_count; ++i) {
        double value = 0.0;
        for (std::size_t grid = 0; grid < grid_points_bohr.size(); ++grid) {
            const double rp = norm3(subtract3(atom_coordinates_bohr[static_cast<std::size_t>(i)], grid_points_bohr[grid]));
            value += mep[grid] / rp;
        }
        matrix_b[static_cast<std::size_t>(i)] = value;
    }
    matrix_b.back() = static_cast<double>(charge);
    result.timings["assembly"] = std::chrono::duration<double>(std::chrono::steady_clock::now() - assembly_start).count();

    const auto stage1_start = std::chrono::steady_clock::now();
    auto q_solution = solve_linear_system(matrix_a, matrix_b);
    std::vector<double> q(q_solution.begin(), q_solution.end() - 1);
    result.esp_charges = force_equivalence_q(q, extra_equivalence);
    if (only_esp) {
        result.stage1_charges = result.esp_charges;
        result.final_charges = result.esp_charges;
        result.timings["stage1"] = std::chrono::duration<double>(std::chrono::steady_clock::now() - stage1_start).count();
        result.timings["stage2"] = 0.0;
        result.timings["total"] = std::chrono::duration<double>(std::chrono::steady_clock::now() - total_start).count();
        return result;
    }

    q = resp_scf_kernel(assign, atom_count, a1, 0.1, matrix_a, matrix_a0, matrix_b, q);
    result.stage1_charges = force_equivalence_q(q, extra_equivalence);
    result.timings["stage1"] = std::chrono::duration<double>(std::chrono::steady_clock::now() - stage1_start).count();
    if (!two_stage) {
        result.final_charges = result.stage1_charges;
        result.timings["stage2"] = 0.0;
        result.timings["total"] = std::chrono::duration<double>(std::chrono::steady_clock::now() - total_start).count();
        return result;
    }

    const auto stage2_start = std::chrono::steady_clock::now();
    std::vector<std::vector<int>> tofit_second;
    std::vector<int> fit_group;
    int sublength = 0;
    find_tofit_second(assign, tofit_second, fit_group, sublength);
    correct_extra_equivalence(tofit_second, fit_group, sublength, extra_equivalence);
    if (!tofit_second.empty()) {
        const int total_length = atom_count - sublength + 1 + atom_count - sublength - static_cast<int>(tofit_second.size());
        std::vector<std::vector<double>> a20;
        std::vector<double> b20;
        get_a20_and_b20(total_length, tofit_second, fit_group, sublength, atom_count, matrix_a0, matrix_b, charge, q, a20, b20);
        result.stage2_restrained_groups = restrained_second_stage_groups(assign, tofit_second);
        auto matrix_stage2 = a20;
        auto q_temp_solution = solve_linear_system(matrix_stage2, b20);
        std::vector<double> q_temp(q_temp_solution.begin(), q_temp_solution.end() - 1);
        std::vector<double> q_last_step = q_temp;
        int step = 0;
        double max_delta = 0.0;
        do {
            ++step;
            q_last_step = q_temp;
            for (int i : result.stage2_restrained_groups) {
                matrix_stage2[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] =
                    a20[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] +
                    a2 / std::sqrt(q_last_step[static_cast<std::size_t>(i)] * q_last_step[static_cast<std::size_t>(i)] + 0.1 * 0.1);
            }
            auto solved = solve_linear_system(matrix_stage2, b20);
            q_temp.assign(solved.begin(), solved.end() - 1);
            max_delta = 0.0;
            for (std::size_t i = 0; i < q_temp.size(); ++i) {
                max_delta = std::max(max_delta, std::abs(q_temp[i] - q_last_step[i]));
            }
        } while (step == 1 || max_delta > 1e-4);
        for (std::size_t i = 0; i < tofit_second.size(); ++i) {
            for (int atom : tofit_second[i]) {
                q[static_cast<std::size_t>(atom)] = q_temp[i];
            }
        }
    }
    result.final_charges = force_equivalence_q(q, extra_equivalence);
    result.timings["stage2"] = std::chrono::duration<double>(std::chrono::steady_clock::now() - stage2_start).count();
    result.timings["total"] = std::chrono::duration<double>(std::chrono::steady_clock::now() - total_start).count();
    return result;
}

std::vector<double> fit_resp_from_esp_cpp(
    const Assign& assign,
    const std::vector<std::array<double, 3>>& atom_coordinates_bohr,
    const std::vector<double>& nuclear_charges,
    const std::vector<std::array<double, 3>>& grid_points_bohr,
    const std::vector<double>& esp_values_au,
    int charge,
    const std::vector<std::vector<int>>& extra_equivalence,
    double a1,
    double a2,
    bool two_stage,
    bool only_esp
) {
    return fit_resp_from_esp_cpp_debug(
        assign,
        atom_coordinates_bohr,
        nuclear_charges,
        grid_points_bohr,
        esp_values_au,
        charge,
        extra_equivalence,
        a1,
        a2,
        two_stage,
        only_esp
    ).final_charges;
}

}  // namespace xpongecpp
