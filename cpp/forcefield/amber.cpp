#include "core.hpp"

#include <array>
#include <fstream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace xpongecpp {
namespace {

std::unordered_map<std::string, ResidueType>& templates() {
    static std::unordered_map<std::string, ResidueType> registry;
    return registry;
}

std::mutex& registry_mutex() {
    static std::mutex mutex;
    return mutex;
}

void put_template(ResidueType residue_type) {
    templates().insert_or_assign(residue_type.name(), std::move(residue_type));
}

void add_minimal_protein_template(const std::string& name) {
    ResidueType residue(name);
    residue.add_atom("N", "N", 0.0, 0.0, 0.0, -0.30, 14.01);
    residue.add_atom("CA", "CT", 1.45, 0.0, 0.0, 0.10, 12.01);
    residue.add_atom("C", "C", 2.05, 1.35, 0.0, 0.50, 12.01);
    residue.add_atom("O", "O", 1.45, 2.35, 0.0, -0.50, 16.00);
    residue.add_connectivity("N", "CA");
    residue.add_connectivity("CA", "C");
    residue.add_connectivity("C", "O");
    put_template(std::move(residue));
}

}  // namespace

void register_ff14sb() {
    std::scoped_lock lock(registry_mutex());
    constexpr std::array<const char*, 30> names{
        "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN",
        "GLU", "GLY", "HID", "HIE", "HIP", "HIS", "ILE", "LEU", "LYS", "MET",
        "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ACE", "NME", "OXT",
    };
    for (const auto* name : names) {
        if (templates().find(name) == templates().end()) {
            add_minimal_protein_template(name);
        }
    }
}

void register_tip3p() {
    std::scoped_lock lock(registry_mutex());

    ResidueType wat("WAT");
    wat.add_atom("O", "OW", 0.0000, 0.0000, 0.0000, -0.834, 16.0);
    wat.add_atom("H1", "HW", 0.9572, 0.0000, 0.0000, 0.417, 1.008);
    wat.add_atom("H2", "HW", -0.2390, 0.9270, 0.0000, 0.417, 1.008);
    wat.add_connectivity("O", "H1");
    wat.add_connectivity("O", "H2");
    put_template(std::move(wat));

    ResidueType na("NA");
    na.add_atom("NA", "Na+", 0.0, 0.0, 0.0, 1.0, 22.99);
    put_template(std::move(na));

    ResidueType cl("CL");
    cl.add_atom("CL", "Cl-", 0.0, 0.0, 0.0, -1.0, 35.45);
    put_template(std::move(cl));
}

void register_residue_templates_from_mol2_text(const std::string& text) {
    Molecule molecule = load_mol2_text(text);
    std::scoped_lock lock(registry_mutex());
    for (const auto& residue : molecule.residues) {
        ResidueType residue_type(residue.name);
        for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
            const auto& atom = molecule.atoms[residue.atom_begin + local];
            residue_type.add_atom(atom.name, atom.type, atom.x, atom.y, atom.z, atom.charge, atom.mass);
        }
        put_template(std::move(residue_type));
    }
}

void register_residue_templates_from_mol2_file(const std::filesystem::path& filename) {
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open mol2 template file: " + filename.string());
    }
    std::ostringstream buffer;
    buffer << input.rdbuf();
    register_residue_templates_from_mol2_text(buffer.str());
}

bool has_template(const std::string& name) {
    std::scoped_lock lock(registry_mutex());
    return templates().find(name) != templates().end();
}

std::size_t template_atom_count(const std::string& name) {
    std::scoped_lock lock(registry_mutex());
    const auto it = templates().find(name);
    if (it == templates().end()) {
        throw std::out_of_range("residue template not found: " + name);
    }
    return it->second.atom_count();
}

const ResidueType& get_residue_template(const std::string& name) {
    std::scoped_lock lock(registry_mutex());
    const auto it = templates().find(name);
    if (it == templates().end()) {
        throw std::out_of_range("residue template not found: " + name);
    }
    return it->second;
}

Molecule get_template_molecule(const std::string& name) {
    const auto& residue_type = get_residue_template(name);
    Molecule molecule(name);
    molecule.append_residue_from_type(residue_type, 0.0, 0.0, 0.0);
    if (!molecule.atoms.empty()) {
        molecule.set_box_padding(0.0, false);
    }
    return molecule;
}

}  // namespace xpongecpp
