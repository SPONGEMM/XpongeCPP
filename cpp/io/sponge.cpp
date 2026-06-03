#include "core.hpp"
#include "sponge_writers.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

namespace xpongecpp {
namespace {

class IndexDisjointSet {
public:
    explicit IndexDisjointSet(std::size_t count) : parent_(count), size_(count, 1) {
        for (std::size_t index = 0; index < count; ++index) {
            parent_[index] = index;
        }
    }

    std::size_t find(std::size_t index) {
        if (parent_[index] != index) {
            parent_[index] = find(parent_[index]);
        }
        return parent_[index];
    }

    void unite(std::size_t lhs, std::size_t rhs) {
        lhs = find(lhs);
        rhs = find(rhs);
        if (lhs == rhs) {
            return;
        }
        if (size_[lhs] < size_[rhs]) {
            std::swap(lhs, rhs);
        }
        parent_[rhs] = lhs;
        size_[lhs] += size_[rhs];
    }

private:
    std::vector<std::size_t> parent_;
    std::vector<std::size_t> size_;
};

std::filesystem::path output_path(const std::filesystem::path& dirname, const std::string& prefix,
                                  const std::string& key) {
    return dirname / (prefix + "_" + key + ".txt");
}

void remember(std::unordered_map<std::string, std::filesystem::path>& outputs, const std::string& key,
              const std::filesystem::path& path) {
    outputs[key] = path;
}

struct OutputBuffer {
    std::string key;
    std::string text;
};

void write_output_buffer(const std::filesystem::path& dirname, const std::string& prefix,
                         std::unordered_map<std::string, std::filesystem::path>& outputs,
                         const OutputBuffer& buffer) {
    const auto path = output_path(dirname, prefix, buffer.key);
    std::ofstream out(path);
    const bool keep_trailing_newline =
        buffer.key == "atom_name" || buffer.key == "atom_type_name" ||
        buffer.key == "resname" || buffer.key == "exclude";
    if (!keep_trailing_newline && !buffer.text.empty() && buffer.text.back() == '\n') {
        out << buffer.text.substr(0, buffer.text.size() - 1);
    } else {
        out << buffer.text;
    }
    remember(outputs, buffer.key, path);
}

std::string formatted_scientific(double value) {
    std::ostringstream out;
    out << std::scientific << std::setw(16) << std::setprecision(7) << value;
    return out.str();
}

std::string python_float(double value) {
    std::ostringstream out;
    out << std::defaultfloat << std::setprecision(12) << value;
    std::string text = out.str();
    if (text.find('.') == std::string::npos && text.find('e') == std::string::npos &&
        text.find('E') == std::string::npos) {
        text += ".0";
    }
    return text;
}

const char* ryckaert_bellemans_listed_definition() {
    return R"([[[ Ryckaert_Bellemans ]]]
[[ parameters ]]
int atom_a, int atom_b, int atom_c, int atom_d, float c0, float c1, float c2, float c3, float c4, float c5
[[ potential ]]
SADfloat<15> cphi = cosf(phi_abcd - CONSTANT_Pi);
SADfloat<15> cphi2 = cphi * cphi;
SADfloat<15> cphi3 = cphi2 * cphi;
SADfloat<15> cphi4 = cphi3 * cphi;
SADfloat<15> cphi5 = cphi4 * cphi;
E = c0 + c1 * cphi + c2 * cphi2 + c3 * cphi3 + c4 * cphi4 + c5 * cphi5;
[[ end ]]
)";
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
    if (lj_combining_rule() == LJCombiningRule::GoodHope) {
        const auto rprod = 4.0 * lj1->second * lj2->second;
        const auto r3 = std::pow(rprod, 3.0);
        return {epsilon * r3 * r3, epsilon * 2.0 * r3};
    }
    const auto rsum = lj1->second + lj2->second;
    const auto r6 = std::pow(rsum, 6.0);
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

bool reorder_residues_by_linked_components(Molecule& molecule) {
    if (molecule.residues.size() < 2 || molecule.residue_links.empty()) {
        return false;
    }

    IndexDisjointSet components(molecule.residues.size());
    for (const auto& link : molecule.residue_links) {
        if (link.atom1 >= molecule.atoms.size() || link.atom2 >= molecule.atoms.size()) {
            throw std::invalid_argument("residue link atom index out of range");
        }
        components.unite(molecule.atoms[link.atom1].residue, molecule.atoms[link.atom2].residue);
    }

    std::unordered_map<std::size_t, double> root_to_sort_key;
    root_to_sort_key.reserve(molecule.residues.size());
    std::vector<double> residue_sort_keys(molecule.residues.size(), 0.0);
    double next_sort_key = 0.0;
    bool already_contiguous = true;
    double previous_key = -1.0;
    for (ResidueId residue_id = 0; residue_id < molecule.residues.size(); ++residue_id) {
        const auto root = components.find(residue_id);
        const auto [it, inserted] = root_to_sort_key.emplace(root, next_sort_key);
        if (inserted) {
            next_sort_key += 1.0;
        }
        const double key = it->second;
        residue_sort_keys[residue_id] = key;
        if (key < previous_key) {
            already_contiguous = false;
        }
        previous_key = key;
    }
    if (already_contiguous) {
        return false;
    }

    molecule.replace_residues({}, residue_sort_keys, true);
    return true;
}

void check_sponge_atom_components_are_contiguous(const Molecule& molecule, const Topology& topology) {
    if (molecule.atoms.empty()) {
        return;
    }

    IndexDisjointSet components(molecule.atoms.size());
    for (const auto& bond : topology.bonds) {
        if (bond.k != 0.0) {
            components.unite(bond.atom1, bond.atom2);
        }
    }

    struct ComponentRange {
        AtomId min_atom{std::numeric_limits<AtomId>::max()};
        AtomId max_atom{0};
        std::size_t count{0};
    };
    std::unordered_map<std::size_t, ComponentRange> ranges;
    ranges.reserve(molecule.atoms.size());
    for (AtomId atom_id = 0; atom_id < molecule.atoms.size(); ++atom_id) {
        auto& range = ranges[components.find(atom_id)];
        range.min_atom = std::min(range.min_atom, atom_id);
        range.max_atom = std::max(range.max_atom, atom_id);
        range.count += 1;
    }

    for (const auto& [root, range] : ranges) {
        (void)root;
        if (static_cast<std::size_t>(range.max_atom - range.min_atom + 1) != range.count) {
            throw std::runtime_error(
                "Atoms in the same molecule must be continuous for SPONGE input; "
                "please reorder residues or atoms before export.");
        }
    }
}

}  // namespace

std::unordered_map<std::string, std::filesystem::path> save_sponge_input(Molecule& input_molecule,
                                                                         const std::string& prefix,
                                                                         const std::filesystem::path& dirname) {
    reorder_residues_by_linked_components(input_molecule);
    std::optional<Molecule> molecule_with_generated_cmaps;
    if (input_molecule.cmaps.empty() && has_amber_cmap_parameters()) {
        molecule_with_generated_cmaps = input_molecule;
        apply_amber_cmaps(*molecule_with_generated_cmaps);
    }
    const Molecule& molecule = molecule_with_generated_cmaps ? *molecule_with_generated_cmaps : input_molecule;
    if (!molecule.validate()) {
        throw std::invalid_argument("cannot export invalid molecule");
    }
    const auto topology = build_topology(molecule);
    check_sponge_atom_components_are_contiguous(molecule, topology);
    const std::string actual_prefix = prefix.empty() ? molecule.name : prefix;
    std::filesystem::create_directories(dirname);
    std::unordered_map<std::string, std::filesystem::path> outputs;

    std::vector<std::pair<std::string, std::future<std::optional<OutputBuffer>>>> futures;
    futures.reserve(12);
    futures.emplace_back("residue", std::async(std::launch::async, [&]() -> std::optional<OutputBuffer> {
        std::ostringstream out;
        out << molecule.atoms.size() << " " << molecule.residues.size() << "\n";
        for (const auto& residue : molecule.residues) {
            out << residue.atom_count << "\n";
        }
        return OutputBuffer{"residue", out.str()};
    }));
    futures.emplace_back("resname", std::async(std::launch::async, [&]() -> std::optional<OutputBuffer> {
        std::ostringstream out;
        out << molecule.residues.size() << "\n";
        for (const auto& residue : molecule.residues) {
            out << residue.name << "\n";
        }
        return OutputBuffer{"resname", out.str()};
    }));
    futures.emplace_back("atom_name", std::async(std::launch::async, [&]() -> std::optional<OutputBuffer> {
        std::ostringstream out;
        out << molecule.atoms.size() << "\n";
        for (const auto& atom : molecule.atoms) {
            out << atom.name << "\n";
        }
        return OutputBuffer{"atom_name", out.str()};
    }));
    futures.emplace_back("atom_type_name", std::async(std::launch::async, [&]() -> std::optional<OutputBuffer> {
        std::ostringstream out;
        out << molecule.atoms.size() << "\n";
        for (const auto& atom : molecule.atoms) {
            out << atom.type << "\n";
        }
        return OutputBuffer{"atom_type_name", out.str()};
    }));
    futures.emplace_back("mass", std::async(std::launch::async, [&]() -> std::optional<OutputBuffer> {
        std::ostringstream out;
        out << molecule.atoms.size() << "\n";
        out << std::fixed << std::setprecision(3);
        for (const auto& atom : molecule.atoms) {
            out << atom.mass << "\n";
        }
        return OutputBuffer{"mass", out.str()};
    }));
    futures.emplace_back("charge", std::async(std::launch::async, [&]() -> std::optional<OutputBuffer> {
        std::ostringstream out;
        out << molecule.atoms.size() << "\n";
        out << std::fixed << std::setprecision(6);
        for (const auto& atom : molecule.atoms) {
            out << atom.charge * 18.2223 << "\n";
        }
        return OutputBuffer{"charge", out.str()};
    }));
    futures.emplace_back("coordinate", std::async(std::launch::async, [&]() -> std::optional<OutputBuffer> {
        std::ostringstream out;
        out << molecule.atoms.size() << "\n";
        out << std::fixed << std::setprecision(6);
        const auto shift = sponge_coordinate_shift_for_export(molecule);
        for (const auto& atom : molecule.atoms) {
            out << atom.x + shift[0] << " " << atom.y + shift[1] << " " << atom.z + shift[2] << "\n";
        }
        const auto box = sponge_coordinate_box_for_export(molecule);
        out << box[0] << " " << box[1] << " " << box[2] << " " << box[3] << " " << box[4] << " " << box[5] << "\n";
        return OutputBuffer{"coordinate", out.str()};
    }));
    futures.emplace_back("LJ", std::async(std::launch::async, [&]() -> std::optional<OutputBuffer> {
        if (molecule.write_lj_soft_core) {
            return std::nullopt;
        }
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

        std::ostringstream out;
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
        return OutputBuffer{"LJ", out.str()};
    }));
    futures.emplace_back("bond", std::async(std::launch::async, [&]() -> std::optional<OutputBuffer> {
        struct BondRow {
            AtomId atom1{0};
            AtomId atom2{0};
            double k{0.0};
            double length{0.0};
        };
        std::vector<BondRow> rows;
        rows.reserve(topology.bonds.size());
        for (const auto& bond : topology.bonds) {
            rows.push_back({bond.atom1, bond.atom2, bond.k, bond.length});
        }
        std::sort(rows.begin(), rows.end(), [](const auto& left, const auto& right) {
            return std::tie(left.atom1, left.atom2) < std::tie(right.atom1, right.atom2);
        });
        std::ostringstream out;
        out << rows.size() << "\n";
        out << std::fixed << std::setprecision(6);
        for (const auto& bond : rows) {
            out << bond.atom1 << " " << bond.atom2 << " " << bond.k << " " << bond.length << "\n";
        }
        return OutputBuffer{"bond", out.str()};
    }));
    futures.emplace_back("angle", std::async(std::launch::async, [&]() -> std::optional<OutputBuffer> {
        struct AngleRow {
            AtomId atom1{0};
            AtomId atom2{0};
            AtomId atom3{0};
            double k{0.0};
            double theta{0.0};
        };
        std::vector<AngleRow> rows;
        rows.reserve(topology.angles.size());
        for (const auto& angle : topology.angles) {
            rows.push_back({angle.atom1, angle.atom2, angle.atom3, angle.k, angle.theta});
        }
        std::sort(rows.begin(), rows.end(), [](const auto& left, const auto& right) {
            return std::tie(left.atom1, left.atom2, left.atom3) < std::tie(right.atom1, right.atom2, right.atom3);
        });
        std::ostringstream out;
        out << rows.size() << "\n";
        out << std::fixed << std::setprecision(6);
        for (const auto& angle : rows) {
            out << angle.atom1 << " " << angle.atom2 << " " << angle.atom3 << " " << angle.k << " " << angle.theta
                << "\n";
        }
        return OutputBuffer{"angle", out.str()};
    }));
    futures.emplace_back("dihedral", std::async(std::launch::async, [&]() -> std::optional<OutputBuffer> {
        struct DihedralRow {
            AtomId atom1{0};
            AtomId atom2{0};
            AtomId atom3{0};
            AtomId atom4{0};
            int periodicity{1};
            double k{0.0};
            double phase{0.0};
        };
        std::vector<DihedralRow> rows;
        rows.reserve(topology.dihedrals.size());
        for (const auto& dihedral : topology.dihedrals) {
            rows.push_back({dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4,
                            dihedral.periodicity, dihedral.k, dihedral.phase});
        }
        std::sort(rows.begin(), rows.end(), [](const auto& left, const auto& right) {
            return std::tie(left.atom1, left.atom2, left.atom3, left.atom4, left.periodicity, left.k, left.phase) <
                   std::tie(right.atom1, right.atom2, right.atom3, right.atom4, right.periodicity, right.k, right.phase);
        });
        std::ostringstream out;
        out << rows.size() << "\n";
        out << std::fixed << std::setprecision(6);
        for (const auto& dihedral : rows) {
            out << dihedral.atom1 << " " << dihedral.atom2 << " " << dihedral.atom3 << " " << dihedral.atom4
                << " " << dihedral.periodicity << " " << dihedral.k << " " << dihedral.phase << "\n";
        }
        return OutputBuffer{"dihedral", out.str()};
    }));
    futures.emplace_back("exclude", std::async(std::launch::async, [&]() -> std::optional<OutputBuffer> {
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
        std::ostringstream out;
        out << topology.exclusions.size() << " " << total_exclusions << "\n";
        for (const auto& exclusions : upper_exclusions) {
            out << exclusions.size();
            out << " ";
            for (std::size_t i = 0; i < exclusions.size(); ++i) {
                if (i != 0) {
                    out << " ";
                }
                out << exclusions[i];
            }
            out << "\n";
        }
        return OutputBuffer{"exclude", out.str()};
    }));
    futures.emplace_back("nb14", std::async(std::launch::async, [&]() -> std::optional<OutputBuffer> {
        std::ostringstream out;
        out << topology.nb14s.size() << "\n";
        out << std::fixed << std::setprecision(6);
        for (const auto& nb14 : topology.nb14s) {
            out << nb14.atom1 << " " << nb14.atom2 << " " << nb14.k_lj << " " << nb14.k_ee << "\n";
        }
        return OutputBuffer{"nb14", out.str()};
    }));
    for (auto& future : futures) {
        const auto buffer = future.second.get();
        if (buffer) {
            write_output_buffer(dirname, actual_prefix, outputs, *buffer);
        }
    }
    if (!molecule.virtual_atoms.empty()) {
        const auto path = output_path(dirname, actual_prefix, "virtual_atom");
        std::ofstream out(path);
        out << std::fixed << std::setprecision(6);
        for (const auto& vatom : molecule.virtual_atoms) {
            out << "2 " << vatom.virtual_atom << " " << vatom.atom0 << " " << vatom.atom1 << " "
                << vatom.atom2 << " " << vatom.k1 << " " << vatom.k2 << "\n";
        }
        remember(outputs, "virtual_atom", path);
    }
    if (!molecule.harmonic_impropers.empty()) {
        struct ImproperRow {
            AtomId atom0{0};
            AtomId atom1{0};
            AtomId atom2{0};
            AtomId atom3{0};
            double k{0.0};
            double phi0{0.0};
        };
        std::vector<ImproperRow> rows;
        rows.reserve(molecule.harmonic_impropers.size());
        for (const auto& improper : molecule.harmonic_impropers) {
            if (improper.k != 0.0) {
                rows.push_back({improper.atom2, improper.atom0, improper.atom1, improper.atom3,
                                improper.k, improper.phi0});
            }
        }
        if (!rows.empty()) {
            std::sort(rows.begin(), rows.end(), [](const auto& left, const auto& right) {
                return std::tuple(left.atom0, left.atom1, left.atom2, left.atom3, left.k, left.phi0) <
                       std::tuple(right.atom0, right.atom1, right.atom2, right.atom3, right.k, right.phi0);
            });
            const auto path = output_path(dirname, actual_prefix, "improper_dihedral");
            std::ofstream out(path);
            out << rows.size() << "\n";
            out << std::fixed << std::setprecision(6);
            for (const auto& row : rows) {
                out << row.atom0 << " " << row.atom1 << " " << row.atom2 << " " << row.atom3 << " "
                    << row.k << " " << row.phi0 << "\n";
            }
            remember(outputs, "improper_dihedral", path);
        }
    }
    if (!molecule.nb14_extras.empty()) {
        struct NB14ExtraRow {
            AtomId atom1{0};
            AtomId atom2{0};
            double a{0.0};
            double b{0.0};
            double kee{0.0};
        };
        std::vector<NB14ExtraRow> rows;
        rows.reserve(molecule.nb14_extras.size());
        for (const auto& nb14 : molecule.nb14_extras) {
            if (nb14.a != 0.0 || nb14.b != 0.0 || nb14.kee != 0.0) {
                const auto atom1 = std::min(nb14.atom1, nb14.atom2);
                const auto atom2 = std::max(nb14.atom1, nb14.atom2);
                rows.push_back({atom1, atom2, nb14.a, nb14.b, nb14.kee});
            }
        }
        if (!rows.empty()) {
            std::sort(rows.begin(), rows.end(), [](const auto& left, const auto& right) {
                return std::tie(left.atom1, left.atom2) < std::tie(right.atom1, right.atom2);
            });
            const auto path = output_path(dirname, actual_prefix, "nb14_extra");
            std::ofstream out(path);
            out << rows.size() << "\n";
            out << std::scientific << std::setprecision(6);
            for (const auto& row : rows) {
                out << row.atom1 << " " << row.atom2 << " " << row.a << " " << row.b << " " << row.kee << "\n";
            }
            remember(outputs, "nb14_extra", path);
        }
    }
    if (!molecule.cmaps.empty()) {
        struct CMapRow {
            AtomId atom0{0};
            AtomId atom1{0};
            AtomId atom2{0};
            AtomId atom3{0};
            AtomId atom4{0};
            std::uint32_t type{0};
        };
        std::vector<CMapRow> rows;
        rows.reserve(molecule.cmaps.size());
        std::vector<std::uint32_t> used_types;
        std::unordered_map<std::uint32_t, std::uint32_t> output_type_by_source_type;
        for (const auto& cmap : molecule.cmaps) {
            auto it = output_type_by_source_type.find(cmap.type);
            if (it == output_type_by_source_type.end()) {
                const auto output_type = static_cast<std::uint32_t>(used_types.size());
                it = output_type_by_source_type.emplace(cmap.type, output_type).first;
                used_types.push_back(cmap.type);
            }
            rows.push_back({cmap.atom0, cmap.atom1, cmap.atom2, cmap.atom3, cmap.atom4, it->second});
        }
        std::sort(rows.begin(), rows.end(), [](const auto& left, const auto& right) {
            return std::tuple(left.atom0, left.atom1, left.atom2, left.atom3, left.atom4) <
                   std::tuple(right.atom0, right.atom1, right.atom2, right.atom3, right.atom4);
        });
        const auto path = output_path(dirname, actual_prefix, "cmap");
        std::ofstream out(path);
        out << molecule.cmaps.size() << " " << used_types.size() << "\n";
        for (const auto type_id : used_types) {
            const auto& type = molecule.cmap_types[type_id];
            out << type.resolution << " ";
        }
        out << "\n";
        out << std::fixed << std::setprecision(6);
        for (const auto type_id : used_types) {
            const auto& type = molecule.cmap_types[type_id];
            for (std::size_t i = 0; i < type.parameters.size(); ++i) {
                out << type.parameters[i] << " ";
                if ((i + 1) % type.resolution == 0) {
                    out << "\n";
                }
            }
            out << "\n";
        }
        for (const auto& row : rows) {
            out << row.atom0 << " " << row.atom1 << " " << row.atom2 << " " << row.atom3 << " "
                << row.atom4 << " " << row.type << "\n";
        }
        remember(outputs, "cmap", path);
    }
    if (!molecule.urey_bradleys.empty()) {
        struct UreyBradleyRow {
            AtomId atom0{0};
            AtomId atom1{0};
            AtomId atom2{0};
            double k{0.0};
            double b{0.0};
            double k_ub{0.0};
            double r13{0.0};
        };
        std::vector<UreyBradleyRow> rows;
        rows.reserve(molecule.urey_bradleys.size());
        for (const auto& angle : molecule.urey_bradleys) {
            if (angle.k != 0.0 || angle.k_ub != 0.0) {
                if (angle.atom0 > angle.atom2) {
                    rows.push_back({angle.atom2, angle.atom1, angle.atom0, angle.k, angle.b, angle.k_ub, angle.r13});
                } else {
                    rows.push_back({angle.atom0, angle.atom1, angle.atom2, angle.k, angle.b, angle.k_ub, angle.r13});
                }
            }
        }
        if (!rows.empty()) {
            std::sort(rows.begin(), rows.end(), [](const auto& left, const auto& right) {
                return std::tie(left.atom0, left.atom1, left.atom2) <
                       std::tie(right.atom0, right.atom1, right.atom2);
            });
            const auto path = output_path(dirname, actual_prefix, "urey_bradley");
            std::ofstream out(path);
            out << rows.size() << "\n";
            out << std::fixed << std::setprecision(6);
            for (const auto& row : rows) {
                out << row.atom0 << " " << row.atom1 << " " << row.atom2 << " " << row.k << " " << row.b
                    << " " << row.k_ub << " " << row.r13 << "\n";
            }
            remember(outputs, "urey_bradley", path);
        }
    }
    if (!molecule.ryckaert_bellemans.empty()) {
        struct RyckaertBellemansRow {
            AtomId atom0{0};
            AtomId atom1{0};
            AtomId atom2{0};
            AtomId atom3{0};
            std::array<double, 6> coefficients{};
        };
        std::vector<RyckaertBellemansRow> rows;
        rows.reserve(molecule.ryckaert_bellemans.size());
        for (const auto& dihedral : molecule.ryckaert_bellemans) {
            const std::array<double, 6> coefficients{dihedral.c0, dihedral.c1, dihedral.c2,
                                                     dihedral.c3, dihedral.c4, dihedral.c5};
            const bool nonzero = std::any_of(coefficients.begin(), coefficients.end(),
                                             [](double value) { return value != 0.0; });
            if (!nonzero) {
                continue;
            }
            if (dihedral.atom0 > dihedral.atom3) {
                rows.push_back({dihedral.atom3, dihedral.atom2, dihedral.atom1, dihedral.atom0, coefficients});
            } else {
                rows.push_back({dihedral.atom0, dihedral.atom1, dihedral.atom2, dihedral.atom3, coefficients});
            }
        }
        if (!rows.empty()) {
            std::sort(rows.begin(), rows.end(), [](const auto& left, const auto& right) {
                return std::tuple(left.atom0, left.atom1, left.atom2, left.atom3, left.coefficients) <
                       std::tuple(right.atom0, right.atom1, right.atom2, right.atom3, right.coefficients);
            });
            const auto path = output_path(dirname, actual_prefix, "Ryckaert_Bellemans");
            std::ofstream out(path);
            out << rows.size() << "\n";
            out << std::fixed << std::setprecision(6);
            for (const auto& row : rows) {
                out << row.atom0 << " " << row.atom1 << " " << row.atom2 << " " << row.atom3;
                for (const double coefficient : row.coefficients) {
                    out << " " << coefficient;
                }
                out << "\n";
            }
            remember(outputs, "Ryckaert_Bellemans", path);
        }
    }
    if (!molecule.soft_bonds.empty()) {
        struct SoftBondRow {
            AtomId atom1{0};
            AtomId atom2{0};
            double k{0.0};
            double b{0.0};
            int from_a_or_b{0};
        };
        std::vector<SoftBondRow> rows;
        rows.reserve(molecule.soft_bonds.size());
        for (const auto& bond : molecule.soft_bonds) {
            if (bond.k != 0.0) {
                const auto atom1 = std::min(bond.atom1, bond.atom2);
                const auto atom2 = std::max(bond.atom1, bond.atom2);
                rows.push_back({atom1, atom2, bond.k, bond.b, bond.from_a_or_b});
            }
        }
        if (!rows.empty()) {
            std::sort(rows.begin(), rows.end(), [](const auto& left, const auto& right) {
                return std::tie(left.atom1, left.atom2) < std::tie(right.atom1, right.atom2);
            });
            const auto path = output_path(dirname, actual_prefix, "bond_soft");
            std::ofstream out(path);
            out << rows.size() << "\n";
            out << std::fixed << std::setprecision(6);
            for (const auto& row : rows) {
                out << row.atom1 << " " << row.atom2 << " " << row.k << " " << row.b << " "
                    << row.from_a_or_b << "\n";
            }
            remember(outputs, "bond_soft", path);
        }
    }
    if (!molecule.ryckaert_bellemans.empty() || !molecule.listed_force_definitions.empty()) {
        std::vector<std::string> definitions;
        const auto append_unique = [&](const std::string& definition) {
            if (!definition.empty() &&
                std::find(definitions.begin(), definitions.end(), definition) == definitions.end()) {
                definitions.push_back(definition);
            }
        };
        if (!molecule.ryckaert_bellemans.empty()) {
            append_unique(ryckaert_bellemans_listed_definition());
        }
        for (const auto& definition : molecule.listed_force_definitions) {
            append_unique(definition);
        }
        if (!definitions.empty()) {
            const auto path = output_path(dirname, actual_prefix, "listed_forces");
            std::ofstream out(path);
            for (const auto& definition : definitions) {
                out << definition;
            }
            remember(outputs, "listed_forces", path);
        }
    }
    if (molecule.has_gb_parameters) {
        const auto path = output_path(dirname, actual_prefix, "gb");
        std::ofstream out(path);
        out << molecule.atoms.size() << "\n";
        out << std::fixed << std::setprecision(4);
        for (const auto& atom : molecule.atoms) {
            out << atom.gb_radius << " " << atom.gb_scaler << "\n";
        }
        remember(outputs, "gb", path);
    }
    if (molecule.write_min_bonded_parameters) {
        {
            const auto path = output_path(dirname, actual_prefix, "fake_mass");
            std::ofstream out(path);
            out << molecule.atoms.size() << "\n";
            out << std::fixed << std::setprecision(3);
            for (const auto& atom : molecule.atoms) {
                const double frozen = (atom.mass < 3.999 || atom.zero_lj_atom || atom.bad_coordinate) ? 1.0 : 0.0;
                out << frozen << "\n";
            }
            remember(outputs, "fake_mass", path);
        }
        {
            const auto path = output_path(dirname, actual_prefix, "fake_LJ");
            std::ofstream out(path);
            out << molecule.atoms.size() << " " << 1 << "\n\n";
            out << formatted_scientific(0.0) << " \n\n";
            out << formatted_scientific(0.0) << " \n\n";
            for (std::size_t i = 0; i < molecule.atoms.size(); ++i) {
                out << 0 << "\n";
            }
            remember(outputs, "fake_LJ", path);
        }
        {
            const auto path = output_path(dirname, actual_prefix, "fake_charge");
            std::ofstream out(path);
            out << molecule.atoms.size() << "\n";
            out << std::fixed << std::setprecision(6);
            for (std::size_t i = 0; i < molecule.atoms.size(); ++i) {
                out << 0.0 << "\n";
            }
            remember(outputs, "fake_charge", path);
        }
    }
    if (molecule.write_subsys_division) {
        const auto path = output_path(dirname, actual_prefix, "subsys_division");
        std::ofstream out(path);
        out << molecule.atoms.size() << "\n";
        for (const auto& atom : molecule.atoms) {
            out << atom.subsys << "\n";
        }
        remember(outputs, "subsys_division", path);
    }
    if (!molecule.sw_parameters.empty()) {
        std::vector<std::string> sw_types;
        std::unordered_map<std::string, std::uint32_t> sw_type_index;
        std::vector<std::uint32_t> atom_type_lines;
        atom_type_lines.reserve(molecule.atoms.size());
        for (const auto& atom : molecule.atoms) {
            if (atom.sw_type.empty()) {
                throw std::runtime_error("missing SW type for atom: " + atom.name);
            }
            auto it = sw_type_index.find(atom.sw_type);
            if (it == sw_type_index.end()) {
                const auto index = static_cast<std::uint32_t>(sw_types.size());
                it = sw_type_index.emplace(atom.sw_type, index).first;
                sw_types.push_back(atom.sw_type);
            }
            atom_type_lines.push_back(it->second);
        }
        const auto path = output_path(dirname, actual_prefix, "SW");
        std::ofstream out(path);
        out << molecule.atoms.size() << " " << sw_types.size() << "\n";
        out << "# type1 type2 A B epsilon[kcal/mol] p q a gamma sigma[Angstrom] "
               "(This is the first required comment line)\n";
        for (const auto& type1 : sw_types) {
            for (const auto& type2 : sw_types) {
                const auto key = type1 + "-" + type2;
                const auto parameter = molecule.sw_parameters.find(key);
                if (parameter == molecule.sw_parameters.end()) {
                    throw std::runtime_error("missing SW parameter: " + key);
                }
                out << sw_type_index.at(type1) << " " << sw_type_index.at(type2) << " "
                    << python_float(parameter->second.a_big) << " " << python_float(parameter->second.b_big) << " "
                    << python_float(parameter->second.epsilon) << " " << python_float(parameter->second.p) << " "
                    << python_float(parameter->second.q) << " " << python_float(parameter->second.a) << " "
                    << python_float(parameter->second.gamma) << " " << python_float(parameter->second.sigma) << "\n";
            }
        }
        out << "# type1 type2 type3 lambda epsilon[kcal/mol] b "
               "(This is the second required comment line)\n";
        for (const auto& type1 : sw_types) {
            for (const auto& type2 : sw_types) {
                for (const auto& type3 : sw_types) {
                    const auto key = type1 + "-" + type2 + "-" + type3;
                    const auto parameter = molecule.sw_parameters.find(key);
                    if (parameter == molecule.sw_parameters.end()) {
                        throw std::runtime_error("missing SW parameter: " + key);
                    }
                    out << sw_type_index.at(type1) << " " << sw_type_index.at(type2) << " "
                        << sw_type_index.at(type3) << " " << python_float(parameter->second.lambda) << " "
                        << python_float(parameter->second.epsilon) << " " << python_float(parameter->second.b) << "\n";
                }
            }
        }
        out << "# atom type from the zeroth atom (This is the third required comment line)\n";
        for (const auto atom_type : atom_type_lines) {
            out << atom_type << "\n";
        }
        remember(outputs, "SW", path);
    }
    if (!molecule.edip_parameters.empty()) {
        std::vector<std::string> edip_types;
        std::unordered_map<std::string, std::uint32_t> edip_type_index;
        std::vector<std::uint32_t> atom_type_lines;
        atom_type_lines.reserve(molecule.atoms.size());
        for (const auto& atom : molecule.atoms) {
            if (atom.edip_type.empty()) {
                throw std::runtime_error("missing EDIP type for atom: " + atom.name);
            }
            auto it = edip_type_index.find(atom.edip_type);
            if (it == edip_type_index.end()) {
                const auto index = static_cast<std::uint32_t>(edip_types.size());
                it = edip_type_index.emplace(atom.edip_type, index).first;
                edip_types.push_back(atom.edip_type);
            }
            atom_type_lines.push_back(it->second);
        }
        const auto path = output_path(dirname, actual_prefix, "EDIP");
        std::ofstream out(path);
        out << molecule.atoms.size() << " " << edip_types.size() << "\n";
        out << "# type1 type2 alpha c[A] a[A] A[kcal/mol] B[A] rho beta sigma[A] "
               "(This is the first required comment line)\n";
        for (const auto& type1 : edip_types) {
            for (const auto& type2 : edip_types) {
                const auto key = type1 + "-" + type2;
                const auto parameter = molecule.edip_parameters.find(key);
                if (parameter == molecule.edip_parameters.end()) {
                    throw std::runtime_error("missing EDIP parameter: " + key);
                }
                out << edip_type_index.at(type1) << " " << edip_type_index.at(type2) << " "
                    << python_float(parameter->second.alpha) << " " << python_float(parameter->second.c) << " "
                    << python_float(parameter->second.a) << " " << python_float(parameter->second.a_big) << " "
                    << python_float(parameter->second.b_big) << " " << python_float(parameter->second.rho) << " "
                    << python_float(parameter->second.beta) << " " << python_float(parameter->second.sigma) << "\n";
            }
        }
        out << "# type1 type2 type3 eta gamma[A] l[kcal/mol] Q0 mu u1 u2 u3 u4 "
               "(This is the second required comment line)\n";
        for (const auto& type1 : edip_types) {
            for (const auto& type2 : edip_types) {
                for (const auto& type3 : edip_types) {
                    const auto key = type1 + "-" + type2 + "-" + type3;
                    const auto parameter = molecule.edip_parameters.find(key);
                    if (parameter == molecule.edip_parameters.end()) {
                        throw std::runtime_error("missing EDIP parameter: " + key);
                    }
                    out << edip_type_index.at(type1) << " " << edip_type_index.at(type2) << " "
                        << edip_type_index.at(type3) << " " << python_float(parameter->second.eta) << " "
                        << python_float(parameter->second.gamma) << " " << python_float(parameter->second.lambda)
                        << " " << python_float(parameter->second.q0) << " " << python_float(parameter->second.mu)
                        << " " << python_float(parameter->second.u1) << " " << python_float(parameter->second.u2)
                        << " " << python_float(parameter->second.u3) << " " << python_float(parameter->second.u4)
                        << "\n";
                }
            }
        }
        out << "# atom type from the zeroth atom (This is the third required comment line)\n";
        for (const auto atom_type : atom_type_lines) {
            out << atom_type << "\n";
        }
        remember(outputs, "EDIP", path);
    }
    if (molecule.write_lj_soft_core) {
        std::vector<std::string> lj_types;
        std::unordered_map<std::string, std::uint32_t> lj_type_index;
        std::vector<std::string> lj_type_b;
        std::unordered_map<std::string, std::uint32_t> lj_type_b_index;
        for (const auto& atom : molecule.atoms) {
            const auto lj_type = find_amber_lj_type(atom.type);
            if (lj_type_index.find(lj_type) == lj_type_index.end()) {
                lj_type_index[lj_type] = static_cast<std::uint32_t>(lj_types.size());
                lj_types.push_back(lj_type);
            }
            const auto type_b = atom.lj_type_b.empty() ? atom.type : atom.lj_type_b;
            const auto lj_b = find_amber_lj_type(type_b);
            if (lj_type_b_index.find(lj_b) == lj_type_b_index.end()) {
                lj_type_b_index[lj_b] = static_cast<std::uint32_t>(lj_type_b.size());
                lj_type_b.push_back(lj_b);
            }
        }

        const auto [full_a, full_b] = find_ab_lj(lj_types, true);
        const auto [full_ab, full_bb] = find_ab_lj(lj_type_b, true);
        const auto checks = lj_check_rows(lj_types, full_a, full_b);
        auto same_type = judge_same_lj_type(lj_types, checks);
        const auto real_types = real_lj_types(lj_types, same_type);
        const auto [real_a, real_b] = find_ab_lj(real_types, false);

        const auto checks_b = lj_check_rows(lj_type_b, full_ab, full_bb);
        auto same_type_b = judge_same_lj_type(lj_type_b, checks_b);
        const auto real_types_b = real_lj_types(lj_type_b, same_type_b);
        const auto [real_ab, real_bb] = find_ab_lj(real_types_b, false);

        const auto path = output_path(dirname, actual_prefix, "LJ_soft_core");
        std::ofstream out(path);
        out << molecule.atoms.size() << " " << real_types.size() << " " << real_types_b.size() << "\n\n";
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
        count = 0;
        for (std::size_t i = 0; i < real_types_b.size(); ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
                (void)j;
                out << formatted_scientific(real_ab[count]) << " ";
                ++count;
            }
            out << "\n";
        }
        out << "\n";
        count = 0;
        for (std::size_t i = 0; i < real_types_b.size(); ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
                (void)j;
                out << formatted_scientific(real_bb[count]) << " ";
                ++count;
            }
            out << "\n";
        }
        out << "\n";
        for (const auto& atom : molecule.atoms) {
            const auto lj_type = find_amber_lj_type(atom.type);
            const auto type_b = atom.lj_type_b.empty() ? atom.type : atom.lj_type_b;
            const auto lj_b = find_amber_lj_type(type_b);
            out << same_type[lj_type_index.at(lj_type)] << " "
                << same_type_b[lj_type_b_index.at(lj_b)] << "\n";
        }
        remember(outputs, "LJ_soft_core", path);
    }

    return outputs;
}

}  // namespace xpongecpp
