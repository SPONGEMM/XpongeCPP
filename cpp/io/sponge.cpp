#include "core.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
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

std::pair<double, double> lj_for_type(const std::string& atom_type) {
    if (atom_type == "OW" || atom_type == "O") {
        return {0.1521, 3.1507};
    }
    if (atom_type == "HW" || atom_type == "H") {
        return {0.0, 1.0};
    }
    if (atom_type == "Na+") {
        return {0.1301, 2.35};
    }
    if (atom_type == "Cl-") {
        return {0.1000, 4.40};
    }
    if (atom_type == "C" || atom_type == "CT") {
        return {0.1094, 3.40};
    }
    if (atom_type == "N") {
        return {0.1700, 3.25};
    }
    if (atom_type == "S") {
        return {0.2500, 3.55};
    }
    return {0.1000, 3.00};
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
        out << std::setprecision(10);
        for (const auto& atom : molecule.atoms) {
            out << atom.mass << "\n";
        }
        remember(outputs, "mass", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "charge");
        std::ofstream out(path);
        out << molecule.atoms.size() << "\n";
        out << std::setprecision(10);
        for (const auto& atom : molecule.atoms) {
            out << atom.charge << "\n";
        }
        remember(outputs, "charge", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "coordinate");
        std::ofstream out(path);
        out << molecule.atoms.size() << "\n";
        out << std::fixed << std::setprecision(6);
        for (const auto& atom : molecule.atoms) {
            out << atom.x << " " << atom.y << " " << atom.z << "\n";
        }
        out << molecule.box_length[0] << " " << molecule.box_length[1] << " " << molecule.box_length[2] << " "
            << molecule.box_angle[0] << " " << molecule.box_angle[1] << " " << molecule.box_angle[2] << "\n";
        remember(outputs, "coordinate", path);
    }
    {
        std::map<std::string, std::uint32_t> type_index;
        for (const auto& atom : molecule.atoms) {
            if (type_index.find(atom.type) == type_index.end()) {
                type_index[atom.type] = static_cast<std::uint32_t>(type_index.size());
            }
        }
        const auto path = output_path(dirname, actual_prefix, "LJ");
        std::ofstream out(path);
        out << type_index.size() << "\n";
        for (const auto& [type, index] : type_index) {
            const auto [epsilon, sigma] = lj_for_type(type);
            out << index << " " << type << " " << epsilon << " " << sigma << "\n";
        }
        out << molecule.atoms.size() << "\n";
        for (const auto& atom : molecule.atoms) {
            out << type_index.at(atom.type) << "\n";
        }
        remember(outputs, "LJ", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "bond");
        std::ofstream out(path);
        out << topology.bonds.size() << "\n";
        out << std::setprecision(10);
        for (const auto& bond : topology.bonds) {
            out << bond.atom1 << " " << bond.atom2 << " " << bond.k << " " << bond.length << "\n";
        }
        remember(outputs, "bond", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "angle");
        std::ofstream out(path);
        out << topology.angles.size() << "\n";
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
        for (const auto& dihedral : topology.dihedrals) {
            out << dihedral.atom1 << " " << dihedral.atom2 << " " << dihedral.atom3 << " " << dihedral.atom4
                << " " << dihedral.periodicity << " " << dihedral.k << " " << dihedral.phase << "\n";
        }
        remember(outputs, "dihedral", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "exclude");
        std::ofstream out(path);
        out << topology.exclusions.size() << "\n";
        for (const auto& exclusions : topology.exclusions) {
            out << exclusions.size();
            for (const auto atom_id : exclusions) {
                out << " " << atom_id;
            }
            out << "\n";
        }
        remember(outputs, "exclude", path);
    }
    {
        const auto path = output_path(dirname, actual_prefix, "nb14");
        std::ofstream out(path);
        out << topology.nb14s.size() << "\n";
        for (const auto& nb14 : topology.nb14s) {
            out << nb14.atom1 << " " << nb14.atom2 << " " << nb14.k_lj << " " << nb14.k_ee << "\n";
        }
        remember(outputs, "nb14", path);
    }

    return outputs;
}

}  // namespace xpongecpp
