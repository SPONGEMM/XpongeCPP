#include "core.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace xpongecpp {
namespace {

std::filesystem::path output_path(const std::filesystem::path& dirname, const std::string& prefix,
                                  const std::string& key) {
    return dirname / (prefix + "_" + key + ".txt");
}

void remember(std::unordered_map<std::string, std::filesystem::path>& outputs, const std::string& key,
              const std::filesystem::path& path) {
    outputs[key] = path;
}

std::array<double, 6> coordinate_box_for_export(const Molecule& molecule) {
    std::array<double, 6> box{molecule.box_length[0], molecule.box_length[1], molecule.box_length[2],
                              molecule.box_angle[0], molecule.box_angle[1], molecule.box_angle[2]};
    if (molecule.has_box || molecule.atoms.empty()) {
        return box;
    }

    std::array<double, 3> minv{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::infinity()};
    std::array<double, 3> maxv{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
                               -std::numeric_limits<double>::infinity()};
    for (const auto& atom : molecule.atoms) {
        minv[0] = std::min(minv[0], atom.x);
        minv[1] = std::min(minv[1], atom.y);
        minv[2] = std::min(minv[2], atom.z);
        maxv[0] = std::max(maxv[0], atom.x);
        maxv[1] = std::max(maxv[1], atom.y);
        maxv[2] = std::max(maxv[2], atom.z);
    }
    box[0] = maxv[0] - minv[0] + 6.0;
    box[1] = maxv[1] - minv[1] + 6.0;
    box[2] = maxv[2] - minv[2] + 6.0;
    return box;
}

std::array<double, 3> coordinate_shift_for_export(const Molecule& molecule) {
    if (molecule.has_box || molecule.atoms.empty()) {
        return {0.0, 0.0, 0.0};
    }
    std::array<double, 3> minv{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::infinity()};
    for (const auto& atom : molecule.atoms) {
        minv[0] = std::min(minv[0], atom.x);
        minv[1] = std::min(minv[1], atom.y);
        minv[2] = std::min(minv[2], atom.z);
    }
    return {3.0 - minv[0], 3.0 - minv[1], 3.0 - minv[2]};
}

std::string formatted_scientific(double value) {
    std::ostringstream out;
    out << std::scientific << std::setw(16) << std::setprecision(7) << value;
    return out.str();
}

std::pair<double, double> amber_lj_ab(const std::string& lj_type1, const std::string& lj_type2) {
    const auto lj1 = find_amber_lj_parameter(lj_type1);
    const auto lj2 = find_amber_lj_parameter(lj_type2);
    if (!lj1) {
        throw std::runtime_error("missing Amber LJ parameter for type: " + lj_type1);
    }
    if (!lj2) {
        throw std::runtime_error("missing Amber LJ parameter for type: " + lj_type2);
    }
    const auto epsilon = std::sqrt(lj1->first * lj2->first);
    const auto rmin = lj1->second + lj2->second;
    const auto r6 = std::pow(rmin, 6.0);
    return {epsilon * r6 * r6, epsilon * 2.0 * r6};
}

std::pair<std::vector<double>, std::vector<double>> find_ab_lj(const std::vector<std::string>& lj_types, bool full) {
    std::vector<double> coefficients_a;
    std::vector<double> coefficients_b;
    const std::size_t total = full ? lj_types.size() * lj_types.size() : lj_types.size() * (lj_types.size() + 1) / 2;
    coefficients_a.reserve(total);
    coefficients_b.reserve(total);
    for (std::size_t i = 0; i < lj_types.size(); ++i) {
        const auto j_max = full ? lj_types.size() : i + 1;
        for (std::size_t j = 0; j < j_max; ++j) {
            const auto [a, b] = amber_lj_ab(lj_types[i], lj_types[j]);
            coefficients_a.push_back(a);
            coefficients_b.push_back(b);
        }
    }
    return {coefficients_a, coefficients_b};
}

std::vector<std::string> lj_check_rows(const std::vector<std::string>& lj_types, const std::vector<double>& a,
                                       const std::vector<double>& b) {
    std::vector<std::string> checks(lj_types.size());
    std::size_t count = 0;
    for (std::size_t i = 0; i < lj_types.size(); ++i) {
        std::string row;
        for (std::size_t j = 0; j < lj_types.size(); ++j) {
            (void)j;
            row += formatted_scientific(a[count]) + " ";
            ++count;
        }
        count -= lj_types.size();
        for (std::size_t j = 0; j < lj_types.size(); ++j) {
            (void)j;
            row += formatted_scientific(b[count]) + " ";
            ++count;
        }
        checks[i] = std::move(row);
    }
    return checks;
}

std::vector<std::uint32_t> judge_same_lj_type(const std::vector<std::string>& lj_types,
                                             const std::vector<std::string>& checks) {
    std::vector<std::uint32_t> same_type(lj_types.size());
    for (std::size_t i = 0; i < lj_types.size(); ++i) {
        same_type[i] = static_cast<std::uint32_t>(i);
    }
    for (std::int64_t i = static_cast<std::int64_t>(lj_types.size()) - 1; i >= 0; --i) {
        for (std::size_t j = static_cast<std::size_t>(i) + 1; j < lj_types.size(); ++j) {
            if (checks[static_cast<std::size_t>(i)] == checks[j]) {
                same_type[j] = static_cast<std::uint32_t>(i);
            }
        }
    }
    return same_type;
}

std::vector<std::string> real_lj_types(const std::vector<std::string>& lj_types,
                                       std::vector<std::uint32_t>& same_type) {
    std::vector<std::string> real;
    std::uint32_t to_subtract = 0;
    for (std::size_t i = 0; i < lj_types.size(); ++i) {
        if (same_type[i] == i) {
            real.push_back(lj_types[i]);
            same_type[i] -= to_subtract;
        } else {
            same_type[i] = same_type[same_type[i]];
            ++to_subtract;
        }
    }
    return real;
}

}  // namespace

std::unordered_map<std::string, std::filesystem::path> save_sponge_input(const Molecule& molecule,
                                                                         const std::string& prefix,
                                                                         const std::filesystem::path& dirname) {
    if (!molecule.validate()) {
        throw std::invalid_argument("cannot export invalid molecule");
    }
    const auto topology = build_topology(molecule);
    const std::string actual_prefix = prefix.empty() ? molecule.name : prefix;
    std::filesystem::create_directories(dirname);
    std::unordered_map<std::string, std::filesystem::path> outputs;

    {
        const auto path = output_path(dirname, actual_prefix, "residue");
        std::ofstream out(path);
        out << molecule.atoms.size() << " " << molecule.residues.size() << "\n";
        for (const auto& residue : molecule.residues) {
            out << residue.atom_count << "\n";
        }
        remember(outputs, "residue", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "resname");
        std::ofstream out(path);
        out << molecule.residues.size() << "\n";
        for (const auto& residue : molecule.residues) {
            out << residue.name << "\n";
        }
        remember(outputs, "resname", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "atom_name");
        std::ofstream out(path);
        out << molecule.atoms.size() << "\n";
        for (const auto& atom : molecule.atoms) {
            out << atom.name << "\n";
        }
        remember(outputs, "atom_name", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "atom_type_name");
        std::ofstream out(path);
        out << molecule.atoms.size() << "\n";
        for (const auto& atom : molecule.atoms) {
            out << atom.type << "\n";
        }
        remember(outputs, "atom_type_name", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "mass");
        std::ofstream out(path);
        out << molecule.atoms.size() << "\n";
        out << std::fixed << std::setprecision(3);
        for (const auto& atom : molecule.atoms) {
            out << atom.mass << "\n";
        }
        remember(outputs, "mass", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "charge");
        std::ofstream out(path);
        out << molecule.atoms.size() << "\n";
        out << std::fixed << std::setprecision(6);
        for (const auto& atom : molecule.atoms) {
            out << atom.charge * 18.2223 << "\n";
        }
        remember(outputs, "charge", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "coordinate");
        std::ofstream out(path);
        out << molecule.atoms.size() << "\n";
        out << std::fixed << std::setprecision(6);
        const auto shift = coordinate_shift_for_export(molecule);
        for (const auto& atom : molecule.atoms) {
            out << atom.x + shift[0] << " " << atom.y + shift[1] << " " << atom.z + shift[2] << "\n";
        }
        const auto box = coordinate_box_for_export(molecule);
        out << box[0] << " " << box[1] << " " << box[2] << " " << box[3] << " " << box[4] << " " << box[5] << "\n";
        remember(outputs, "coordinate", path);
    }
    {
        std::vector<std::string> lj_types;
        std::unordered_map<std::string, std::uint32_t> lj_type_index;
        lj_type_index.reserve(32);
        for (const auto& atom : molecule.atoms) {
            const auto lj_type = find_amber_lj_type(atom.type);
            if (lj_type_index.find(lj_type) == lj_type_index.end()) {
                lj_type_index[lj_type] = static_cast<std::uint32_t>(lj_types.size());
                lj_types.push_back(lj_type);
            }
        }
        const auto [full_a, full_b] = find_ab_lj(lj_types, true);
        const auto checks = lj_check_rows(lj_types, full_a, full_b);
        auto same_type = judge_same_lj_type(lj_types, checks);
        const auto real_types = real_lj_types(lj_types, same_type);
        const auto [real_a, real_b] = find_ab_lj(real_types, false);

        const auto path = output_path(dirname, actual_prefix, "LJ");
        std::ofstream out(path);
        out << molecule.atoms.size() << " " << real_types.size() << "\n\n";
        std::size_t count = 0;
        for (std::size_t i = 0; i < real_types.size(); ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
                (void)j;
                out << formatted_scientific(real_a[count]) << " ";
                ++count;
            }
            out << "\n";
        }
        out << "\n";
        count = 0;
        for (std::size_t i = 0; i < real_types.size(); ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
                (void)j;
                out << formatted_scientific(real_b[count]) << " ";
                ++count;
            }
            out << "\n";
        }
        out << "\n";
        for (const auto& atom : molecule.atoms) {
            const auto lj_type = find_amber_lj_type(atom.type);
            out << same_type[lj_type_index.at(lj_type)] << "\n";
        }
        remember(outputs, "LJ", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "bond");
        std::ofstream out(path);
        out << topology.bonds.size() << "\n";
        out << std::fixed << std::setprecision(6);
        for (const auto& bond : topology.bonds) {
            out << bond.atom1 << " " << bond.atom2 << " " << bond.k << " " << bond.length << "\n";
        }
        remember(outputs, "bond", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "angle");
        std::ofstream out(path);
        out << topology.angles.size() << "\n";
        out << std::fixed << std::setprecision(6);
        for (const auto& angle : topology.angles) {
            out << angle.atom1 << " " << angle.atom2 << " " << angle.atom3 << " " << angle.k << " " << angle.theta
                << "\n";
        }
        remember(outputs, "angle", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "dihedral");
        std::ofstream out(path);
        out << topology.dihedrals.size() << "\n";
        out << std::fixed << std::setprecision(6);
        for (const auto& dihedral : topology.dihedrals) {
            out << dihedral.atom1 << " " << dihedral.atom2 << " " << dihedral.atom3 << " " << dihedral.atom4
                << " " << dihedral.periodicity << " " << dihedral.k << " " << dihedral.phase << "\n";
        }
        remember(outputs, "dihedral", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "exclude");
        std::ofstream out(path);
        std::size_t total_exclusions = 0;
        std::vector<std::vector<AtomId>> upper_exclusions(topology.exclusions.size());
        for (AtomId atom_id = 0; atom_id < topology.exclusions.size(); ++atom_id) {
            for (const auto excluded : topology.exclusions[atom_id]) {
                if (excluded > atom_id) {
                    upper_exclusions[atom_id].push_back(excluded);
                }
            }
            std::sort(upper_exclusions[atom_id].begin(), upper_exclusions[atom_id].end());
            total_exclusions += upper_exclusions[atom_id].size();
        }
        out << topology.exclusions.size() << " " << total_exclusions << "\n";
        for (const auto& exclusions : upper_exclusions) {
            out << exclusions.size();
            if (!exclusions.empty()) {
                out << " ";
            }
            for (std::size_t i = 0; i < exclusions.size(); ++i) {
                if (i != 0) {
                    out << " ";
                }
                out << exclusions[i];
            }
            out << "\n";
        }
        remember(outputs, "exclude", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "nb14");
        std::ofstream out(path);
        out << topology.nb14s.size() << "\n";
        out << std::fixed << std::setprecision(6);
        for (const auto& nb14 : topology.nb14s) {
            out << nb14.atom1 << " " << nb14.atom2 << " " << nb14.k_lj << " " << nb14.k_ee << "\n";
        }
        remember(outputs, "nb14", path);
    }

    return outputs;
}

}  // namespace xpongecpp
