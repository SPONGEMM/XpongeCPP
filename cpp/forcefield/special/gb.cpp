#include "core.hpp"

#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace xpongecpp {
namespace {

constexpr AtomId invalid_atom_id = std::numeric_limits<AtomId>::max();

std::pair<double, double> bondi_radius_and_scaler(const std::string& element) {
    if (element == "H") return {1.2, 0.85};
    if (element == "C") return {1.7, 0.72};
    if (element == "N") return {1.55, 0.79};
    if (element == "O") return {1.52, 0.85};
    if (element == "F") return {1.47, 0.88};
    if (element == "P") return {1.8, 0.86};
    if (element == "S") return {1.8, 0.96};
    if (element == "Cl") return {1.75, 0.8};
    if (element == "Br") return {1.85, 0.8};
    if (element == "I") return {1.98, 0.8};
    return {1.5, 0.8};
}

std::pair<double, double> modified_bondi_radius_and_scaler(const std::string& element) {
    if (element == "H") return {1.2, 0.85};
    if (element == "C") return {1.7, 0.72};
    if (element == "N") return {1.55, 0.79};
    if (element == "O") return {1.5, 0.85};
    if (element == "F") return {1.5, 0.88};
    if (element == "Si") return {2.1, 0.8};
    if (element == "P") return {1.85, 0.86};
    if (element == "S") return {1.8, 0.96};
    if (element == "Cl") return {1.7, 0.8};
    if (element == "Br") return {1.85, 0.8};
    if (element == "I") return {1.98, 0.8};
    return {1.5, 0.8};
}

std::string gb_element_from_mass(const Atom& atom) {
    return atom.mass > 0.0 ? guess_element_from_mass(atom.mass) : atom.element;
}

std::string gb_element_from_mass_value(double mass, const std::string& fallback_element = "") {
    return mass > 0.0 ? guess_element_from_mass(mass) : fallback_element;
}

}  // namespace

void Molecule::set_gb_radius(const std::string& radius_set) {
    for (AtomId atom_id = 0; atom_id < atoms.size(); ++atom_id) {
        auto& atom = atoms[atom_id];
        const auto element = gb_element_from_mass(atom);
        if (radius_set == "bondi_radii") {
            const auto [radius, scaler] = bondi_radius_and_scaler(element);
            atom.gb_radius = radius;
            atom.gb_scaler = scaler;
        } else if (radius_set == "modified_bondi_radii") {
            auto [radius, scaler] = modified_bondi_radius_and_scaler(element);
            if (element == "H") {
                AtomId bonded_atom_id = invalid_atom_id;
                for (const auto& bond : explicit_bonds) {
                    if (bond.atom1 == atom_id) {
                        bonded_atom_id = bond.atom2;
                        break;
                    }
                    if (bond.atom2 == atom_id) {
                        bonded_atom_id = bond.atom1;
                        break;
                    }
                }
                const auto residue_id = atom.residue;
                if (bonded_atom_id != invalid_atom_id && bonded_atom_id < atoms.size()) {
                    const auto bonded_element = gb_element_from_mass(atoms[bonded_atom_id]);
                    if (bonded_element == "C" || bonded_element == "N") {
                        radius = 1.3;
                    } else if (bonded_element == "S" || bonded_element == "O" || bonded_element == "H") {
                        radius = 0.8;
                    }
                } else if (residue_id < residues.size()) {
                    const auto& residue = residues[residue_id];
                    if (has_template(residue.name)) {
                        const auto& residue_type = get_residue_template(residue.name);
                        auto template_atom_index = invalid_atom_id;
                        try {
                            template_atom_index = residue_type.atom_index(atom.name);
                        } catch (const std::exception&) {
                            template_atom_index = invalid_atom_id;
                        }
                        if (template_atom_index != invalid_atom_id) {
                            for (const auto& bond : residue_type.bonds()) {
                                AtomId bonded_template = invalid_atom_id;
                                if (bond.atom1 == template_atom_index) {
                                    bonded_template = bond.atom2;
                                } else if (bond.atom2 == template_atom_index) {
                                    bonded_template = bond.atom1;
                                }
                                if (bonded_template == invalid_atom_id) {
                                    continue;
                                }
                                const auto& bonded_atom = residue_type.atoms()[bonded_template];
                                const auto bonded_element =
                                    !bonded_atom.element.empty() ? bonded_atom.element :
                                                                  gb_element_from_mass_value(bonded_atom.mass);
                                if (bonded_element == "C" || bonded_element == "N") {
                                    radius = 1.3;
                                } else if (bonded_element == "S" || bonded_element == "O" || bonded_element == "H") {
                                    radius = 0.8;
                                }
                                break;
                            }
                        }
                    }
                }
            }
            atom.gb_radius = radius;
            atom.gb_scaler = scaler;
        } else {
            throw std::invalid_argument("unknown GB radius set: " + radius_set);
        }
    }
    box_length = {999.0, 999.0, 999.0};
    has_box = true;
    has_gb_parameters = true;
}

}  // namespace xpongecpp
