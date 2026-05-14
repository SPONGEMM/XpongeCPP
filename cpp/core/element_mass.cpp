#include "core.hpp"

#include <array>
#include <cctype>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace xpongecpp {
namespace {

constexpr std::array<const char*, 112> kElements{
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc",
    "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",
    "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
    "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg",
};

constexpr std::array<double, 112> kElementMasses{
    0.00000, 1.00794, 4.00260, 6.941, 9.012182, 10.811,
    12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
    22.989770, 24.3050, 26.981538, 28.0855, 30.973761,
    32.065, 35.453, 39.948, 39.0983, 40.078, 44.955910,
    47.867, 50.9415, 51.9961, 54.938049, 55.845, 58.9332,
    58.6934, 63.546, 65.409, 69.723, 72.64, 74.92160,
    78.96, 79.904, 83.798, 85.4678, 87.62, 88.90585,
    91.224, 92.90638, 95.94, 98.0, 101.07, 102.90550,
    106.42, 107.8682, 112.411, 114.818, 118.710, 121.760,
    127.60, 126.90447, 131.293, 132.90545, 137.327,
    138.9055, 140.116, 140.90765, 144.24, 145.0, 150.36,
    151.964, 157.25, 158.92534, 162.500, 164.93032,
    167.259, 168.93421, 173.04, 174.967, 178.49, 180.9479,
    183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,
    200.59, 204.3833, 207.2, 208.98038, 209.0, 210.0, 222.0,
    223.0, 226.0, 227.0, 232.0381, 231.03588, 238.02891,
    237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0,
    258.0, 259.0, 262.0, 261.0, 262.0, 266.0, 264.0, 269.0,
    268.0, 271.0, 272.0,
};

std::string trim_copy(const std::string& input) {
    const auto first = input.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}

std::string fallback_element_from_name(const std::string& atom_name) {
    std::string letters;
    for (const char c : atom_name) {
        if (std::isalpha(static_cast<unsigned char>(c))) {
            letters.push_back(c);
        }
    }
    if (letters.empty()) {
        return "X";
    }
    std::string element;
    element.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(letters[0]))));
    if (letters.size() > 1 && std::islower(static_cast<unsigned char>(letters[1]))) {
        element.push_back(letters[1]);
    }
    return element;
}

std::string canonicalize_explicit_element(const std::string& explicit_element) {
    const auto trimmed = trim_copy(explicit_element);
    if (trimmed.empty()) {
        return "";
    }
    std::string letters;
    for (const char c : trimmed) {
        if (std::isalpha(static_cast<unsigned char>(c))) {
            letters.push_back(c);
        }
    }
    if (letters.empty()) {
        return "";
    }
    std::string element;
    element.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(letters[0]))));
    if (letters.size() > 1) {
        element.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(letters[1]))));
    }
    return element;
}

std::optional<double> registered_mass_for_type(const std::string& atom_type) {
    if (atom_type.empty()) {
        return std::nullopt;
    }
    if (const auto amber_mass = find_amber_atom_type_mass(atom_type)) {
        return amber_mass;
    }
    if (const auto external_mass = find_external_atom_type_mass(atom_type)) {
        return external_mass;
    }
    return std::nullopt;
}

double resolve_mass_for_atom(const std::string& atom_name, const std::string& atom_type,
                             const std::string& explicit_element, double explicit_mass) {
    if (explicit_mass > 0.0) {
        return explicit_mass;
    }
    if (const auto registered_mass = registered_mass_for_type(atom_type)) {
        return *registered_mass;
    }
    const auto explicit_canonical = canonicalize_explicit_element(explicit_element);
    if (!explicit_canonical.empty()) {
        const auto default_mass = default_mass_for_element(explicit_canonical);
        if (default_mass > 0.0) {
            return default_mass;
        }
    }
    const auto fallback_element = fallback_element_from_name(atom_name);
    return default_mass_for_element(fallback_element);
}

std::string resolve_element_for_atom(const std::string& atom_name, const std::string& explicit_element,
                                     double resolved_mass) {
    const auto explicit_canonical = canonicalize_explicit_element(explicit_element);
    if (!explicit_canonical.empty()) {
        return explicit_canonical;
    }
    if (resolved_mass > 0.0) {
        return guess_element_from_mass(resolved_mass);
    }
    return fallback_element_from_name(atom_name);
}

}  // namespace

double default_mass_for_element(const std::string& element) {
    if (element == "H") return 1.008;
    if (element == "C") return 12.010;
    if (element == "N") return 14.010;
    if (element == "O") return 16.000;
    if (element == "P") return 30.970;
    if (element == "S") return 32.060;
    if (element == "Cl" || element == "CL") return 35.450;
    if (element == "Na" || element == "NA") return 22.990;
    if (element == "K") return 39.100;
    return 0.0;
}

std::string guess_element_from_mass(double mass) {
    std::size_t index = 0;
    if (mass > 0.0 && mass < 3.8) {
        index = 1;
    } else if (mass > 207.85 && mass < 208.99) {
        index = 83;
    } else if (mass > 56.50 && mass < 58.8133) {
        index = 27;
    } else {
        for (std::size_t j = 0; j < kElementMasses.size(); ++j) {
            if (std::abs(mass - kElementMasses[j]) < 0.65) {
                index = j;
                break;
            }
        }
    }
    return kElements[index];
}

std::string guess_element(const std::string& atom_name, const std::string& explicit_element) {
    const auto explicit_canonical = canonicalize_explicit_element(explicit_element);
    if (!explicit_canonical.empty()) {
        return explicit_canonical;
    }
    return fallback_element_from_name(atom_name);
}

void ResidueType::add_atom(const std::string& name, const std::string& type, double x, double y, double z,
                           double charge, double mass) {
    if (atom_name_to_index_.count(name) != 0) {
        throw std::invalid_argument("duplicate atom name in ResidueType: " + name);
    }
    ResidueTypeAtom atom;
    atom.name = name;
    atom.type = type;
    atom.mass = resolve_mass_for_atom(name, type, "", mass);
    atom.element = resolve_element_for_atom(name, "", atom.mass);
    atom.x = x;
    atom.y = y;
    atom.z = z;
    atom.charge = charge;
    atom_name_to_index_[name] = static_cast<std::uint32_t>(atoms_.size());
    atoms_.push_back(std::move(atom));
    ++version_;
}

}  // namespace xpongecpp
