#include "amber_internal.hpp"

#include <array>
#include <mutex>
#include <shared_mutex>

namespace xpongecpp {

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

void register_ff14sb() {
    constexpr std::array<const char*, 30> names{
        "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN",
        "GLU", "GLY", "HID", "HIE", "HIP", "HIS", "ILE", "LEU", "LYS", "MET",
        "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ACE", "NME", "OXT",
    };
    for (const auto* name : names) {
        std::shared_lock lock(registry_mutex());
        const bool exists = templates().find(name) != templates().end();
        lock.unlock();
        if (!exists) {
            add_minimal_protein_template(name);
        }
    }
}

void register_tip3p() {
    ResidueType wat("WAT");
    wat.add_atom("O", "OW", 0.0000, 0.0000, 0.0000, -0.834, 16.0);
    wat.add_atom("H1", "HW", 0.9572, 0.0000, 0.0000, 0.417, 1.008);
    wat.add_atom("H2", "HW", -0.239988, 0.926627, 0.0000, 0.417, 1.008);
    wat.add_connectivity("O", "H1");
    wat.add_connectivity("O", "H2");
    wat.add_connectivity("H1", "H2");
    put_template(std::move(wat));

    ResidueType na("NA");
    na.add_atom("NA", "Na+", 0.0, 0.0, 0.0, 1.0, 22.99);
    put_template(std::move(na));

    ResidueType cl("CL");
    cl.add_atom("CL", "Cl-", 0.0, 0.0, 0.0, -1.0, 35.45);
    put_template(std::move(cl));

    std::unique_lock lock(registry_mutex());
    upsert_lj_atom_type("OW", "OW");
    upsert_lj_atom_type("HW", "HW");
    upsert_lj_parameter("OW", 0.152, 1.7683);
    upsert_lj_parameter("HW", 0.0, 0.0);
    upsert_lj_atom_type("Na+", "Na+");
    upsert_lj_atom_type("Cl-", "Cl-");
}

}  // namespace xpongecpp
