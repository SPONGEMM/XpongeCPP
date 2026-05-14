#include "amber_internal.hpp"

#include <string>
#include <unordered_set>

namespace xpongecpp {

bool is_standard_protein_residue(const std::string& name) {
    static const std::unordered_set<std::string> residues{
        "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN", "GLU", "GLY", "HID", "HIE",
        "HIP", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
    return residues.count(name) != 0;
}

void configure_xponge_residue_links(ResidueType& residue_type) {
    const auto name = residue_type.name();
    if (is_standard_protein_residue(name)) {
        residue_type.set_head("N", 1.3, "CA");
        residue_type.set_tail("C", 1.3, "CA");
        register_pdb_residue_name_mapping("head", name, "N" + name);
        register_pdb_residue_name_mapping("tail", name, "C" + name);
    } else if (name.size() > 1 && name[0] == 'N' && is_standard_protein_residue(name.substr(1))) {
        residue_type.set_tail("C", 1.3, "CA");
    } else if (name.size() > 1 && name[0] == 'C' && is_standard_protein_residue(name.substr(1))) {
        residue_type.set_head("N", 1.3, "CA");
    } else if (name == "ACE") {
        residue_type.set_tail("C", 1.3, "CH3");
    } else if (name == "NME") {
        residue_type.set_head("N", 1.3, "CH3");
    }
    if (name == "CYX" || name == "NCYX" || name == "CCYX") {
        residue_type.set_connect_atom("ssbond", "SG");
    }
    if (name == "HIS") {
        register_his_mapping("HIS", "HID", "HIE", "HIP");
    } else if (name == "NHIS") {
        register_his_mapping("NHIS", "NHID", "NHIE", "NHIP");
    } else if (name == "CHIS") {
        register_his_mapping("CHIS", "CHID", "CHIE", "CHIP");
    }
}

}  // namespace xpongecpp
