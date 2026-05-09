#include "core.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace xpongecpp {
namespace {

std::string trim_copy(const std::string& input) {
    const auto first = input.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}

std::string upper_copy(std::string value) {
    for (auto& ch : value) {
        ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    }
    return value;
}

std::string pdb_string(const std::string& line, std::size_t pos, std::size_t len) {
    if (line.size() <= pos) {
        return "";
    }
    return trim_copy(line.substr(pos, std::min(len, line.size() - pos)));
}

double pdb_float(const std::string& line, std::size_t pos, std::size_t len, const char* field_name) {
    if (line.size() <= pos) {
        throw std::invalid_argument(std::string("missing PDB coordinate field: ") + field_name);
    }
    const auto field = trim_copy(line.substr(pos, std::min(len, line.size() - pos)));
    if (field.empty()) {
        throw std::invalid_argument(std::string("empty PDB coordinate field: ") + field_name);
    }
    try {
        return std::stod(field);
    } catch (const std::exception&) {
        throw std::invalid_argument(std::string("invalid PDB coordinate field ") + field_name + ": " + field);
    }
}

void set_atom_defaults(Atom& atom, const std::string& residue_name) {
    const auto residue_upper = upper_copy(residue_name);
    const auto atom_upper = upper_copy(atom.name);
    if (residue_upper == "WAT" || residue_upper == "HOH" || residue_upper == "TIP3") {
        atom.type = atom.element == "H" ? "HW" : "OW";
        atom.charge = atom.element == "H" ? 0.417 : -0.834;
    } else if (residue_upper == "NA" || residue_upper == "Na+") {
        atom.type = "Na+";
        atom.element = "Na";
        atom.charge = 1.0;
    } else if (residue_upper == "CL" || residue_upper == "Cl-") {
        atom.type = "Cl-";
        atom.element = "Cl";
        atom.charge = -1.0;
    } else if (atom_upper == "OXT") {
        atom.type = "O";
    } else {
        atom.type = atom.element;
    }
    atom.mass = default_mass_for_element(atom.element);
}

std::string pdb_atom_name(const std::string& name) {
    if (name.size() >= 4) {
        return name.substr(0, 4);
    }
    std::ostringstream out;
    out << std::left << std::setw(4) << name;
    return out.str();
}

bool is_terminal_mappable(const Residue& residue) {
    return residue.name != "WAT" && residue.name != "HOH" && residue.name != "NA" && residue.name != "CL" &&
           residue.name.rfind("N", 0) != 0 && residue.name.rfind("C", 0) != 0;
}

void apply_template_atom_properties(Molecule& molecule) {
    for (auto& residue : molecule.residues) {
        if (!has_template(residue.name)) {
            continue;
        }
        const auto& residue_type = get_residue_template(residue.name);
        for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
            auto& atom = molecule.atoms[residue.atom_begin + local];
            try {
                const auto template_index = residue_type.atom_index(atom.name);
                const auto& template_atom = residue_type.atoms()[template_index];
                atom.type = template_atom.type;
                atom.element = template_atom.element;
                atom.charge = template_atom.charge;
                atom.mass = template_atom.mass;
            } catch (const std::exception&) {
            }
        }
    }
}

}  // namespace

Molecule load_pdb_text(const std::string& text) {
    Molecule molecule("PDB");
    std::istringstream input(text);
    std::string line;
    std::string current_residue_key;
    ResidueId current_residue_id = 0;

    while (std::getline(input, line)) {
        if (line.rfind("CRYST1", 0) == 0) {
            try {
                molecule.box_length = {
                    std::stod(trim_copy(line.substr(6, 9))),
                    std::stod(trim_copy(line.substr(15, 9))),
                    std::stod(trim_copy(line.substr(24, 9))),
                };
                molecule.box_angle = {
                    std::stod(trim_copy(line.substr(33, 7))),
                    std::stod(trim_copy(line.substr(40, 7))),
                    std::stod(trim_copy(line.substr(47, 7))),
                };
                molecule.has_box = true;
            } catch (const std::exception&) {
                throw std::invalid_argument("invalid CRYST1 record");
            }
            continue;
        }
        if (!(line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0)) {
            continue;
        }
        const std::string atom_name = pdb_string(line, 12, 4);
        const std::string residue_name = pdb_string(line, 17, 3);
        const std::string chain = pdb_string(line, 21, 1);
        const std::string resseq = pdb_string(line, 22, 4);
        const std::string icode = pdb_string(line, 26, 1);
        const std::string residue_key = chain + ":" + resseq + ":" + icode + ":" + residue_name;

        if (molecule.residues.empty() || residue_key != current_residue_key) {
            current_residue_key = residue_key;
            current_residue_id = static_cast<ResidueId>(molecule.residues.size());
            Residue residue;
            residue.name = residue_name.empty() ? "UNK" : upper_copy(residue_name);
            residue.type_name = residue.name;
            residue.atom_begin = static_cast<AtomId>(molecule.atoms.size());
            residue.atom_count = 0;
            molecule.residues.push_back(residue);
        }

        Atom atom;
        atom.name = atom_name;
        atom.element = guess_element(atom_name, pdb_string(line, 76, 2));
        atom.residue = current_residue_id;
        atom.x = pdb_float(line, 30, 8, "x");
        atom.y = pdb_float(line, 38, 8, "y");
        atom.z = pdb_float(line, 46, 8, "z");
        set_atom_defaults(atom, molecule.residues[current_residue_id].name);
        molecule.atoms.push_back(std::move(atom));
        molecule.residues[current_residue_id].atom_count += 1;
    }
    if (!molecule.residues.empty()) {
        if (is_terminal_mappable(molecule.residues.front())) {
            molecule.residues.front().name = "N" + molecule.residues.front().name;
            molecule.residues.front().type_name = molecule.residues.front().name;
        }
        if (molecule.residues.size() > 1 && is_terminal_mappable(molecule.residues.back())) {
            molecule.residues.back().name = "C" + molecule.residues.back().name;
            molecule.residues.back().type_name = molecule.residues.back().name;
        }
    }
    apply_template_atom_properties(molecule);
    return molecule;
}

void save_pdb(const Molecule& molecule, const std::filesystem::path& filename) {
    if (!molecule.validate()) {
        throw std::invalid_argument("cannot export invalid molecule");
    }
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("failed to open PDB output: " + filename.string());
    }
    out << std::fixed << std::setprecision(3);
    if (molecule.has_box) {
        out << "CRYST1" << std::setw(9) << molecule.box_length[0] << std::setw(9) << molecule.box_length[1]
            << std::setw(9) << molecule.box_length[2] << std::setw(7) << molecule.box_angle[0] << std::setw(7)
            << molecule.box_angle[1] << std::setw(7) << molecule.box_angle[2] << " P 1           1\n";
    }

    std::size_t serial = 1;
    for (std::size_t residue_index = 0; residue_index < molecule.residues.size(); ++residue_index) {
        const auto& residue = molecule.residues[residue_index];
        for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
            const auto& atom = molecule.atoms[residue.atom_begin + local];
            const char* record = (residue.name == "WAT" || residue.name == "NA" || residue.name == "CL") ? "HETATM" : "ATOM  ";
            out << record << std::setw(5) << serial << ' ' << std::left << std::setw(4) << pdb_atom_name(atom.name)
                << std::right << ' ' << std::setw(3) << residue.name << " A" << std::setw(4)
                << ((residue_index % 9999) + 1) << "    " << std::setw(8) << atom.x << std::setw(8) << atom.y
                << std::setw(8) << atom.z << "  1.00  0.00          " << std::setw(2) << atom.element << "\n";
            ++serial;
        }
        out << "TER\n";
    }
    out << "END\n";
}

}  // namespace xpongecpp
