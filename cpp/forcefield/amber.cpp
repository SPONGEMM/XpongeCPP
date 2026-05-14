#include "core.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <mutex>
#include <optional>
#include <shared_mutex>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace xpongecpp {
namespace {

std::unordered_map<std::string, ResidueType>& templates() {
    static std::unordered_map<std::string, ResidueType> registry;
    return registry;
}

std::unordered_map<std::string, Molecule>& molecule_templates() {
    static std::unordered_map<std::string, Molecule> registry;
    return registry;
}

std::shared_mutex& registry_mutex() {
    static std::shared_mutex mutex;
    return mutex;
}

struct ProperParameter {
    std::array<std::string, 4> types;
    std::vector<DihedralTerm> terms;
    std::size_t order{0};
};

struct ImproperParameter {
    std::array<std::string, 4> types;
    DihedralTerm term;
    std::size_t order{0};
};

struct NB14Parameter {
    std::string atom_type1;
    std::string atom_type4;
    NB14Scale scale;
    std::size_t order{0};
};

struct AmberCMapParameter {
    std::uint32_t resolution{0};
    std::vector<double> parameters;
};

struct BondParameter {
    std::string atom_type1;
    std::string atom_type2;
    BondTerm term;
    std::size_t order{0};
};

struct AngleParameter {
    std::array<std::string, 3> types;
    AngleTerm term;
    std::size_t order{0};
};

std::vector<BondParameter>& bond_parameters() {
    static std::vector<BondParameter> parameters;
    return parameters;
}

std::vector<AngleParameter>& angle_parameters() {
    static std::vector<AngleParameter> parameters;
    return parameters;
}

std::vector<ProperParameter>& proper_parameters() {
    static std::vector<ProperParameter> parameters;
    return parameters;
}

std::vector<ImproperParameter>& improper_parameters() {
    static std::vector<ImproperParameter> parameters;
    return parameters;
}

std::vector<NB14Parameter>& nb14_parameters() {
    static std::vector<NB14Parameter> parameters{{"X", "X", {0.5, 0.833333}, 0}};
    return parameters;
}

std::unordered_map<std::string, AmberCMapParameter>& amber_cmap_parameters() {
    static std::unordered_map<std::string, AmberCMapParameter> parameters;
    return parameters;
}

std::unordered_map<std::string, std::string>& lj_type_by_atom_type() {
    static std::unordered_map<std::string, std::string> parameters;
    return parameters;
}

std::unordered_map<std::string, std::pair<double, double>>& lj_parameters() {
    static std::unordered_map<std::string, std::pair<double, double>> parameters;
    return parameters;
}

std::unordered_map<std::string, double>& mass_by_atom_type() {
    static std::unordered_map<std::string, double> parameters;
    return parameters;
}

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

void put_template(ResidueType residue_type) {
    configure_xponge_residue_links(residue_type);
    molecule_templates().erase(residue_type.name());
    templates().insert_or_assign(residue_type.name(), std::move(residue_type));
}

ResidueType residue_type_from_molecule_residue(const Molecule& molecule, const Residue& residue) {
    ResidueType residue_type(residue.name);
    for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
        const auto& atom = molecule.atoms[residue.atom_begin + local];
        residue_type.add_atom(atom.name, atom.type, atom.x, atom.y, atom.z, atom.charge, atom.mass);
    }
    for (const auto& bond : molecule.explicit_bonds) {
        if (bond.atom1 < residue.atom_begin || bond.atom1 >= residue.atom_begin + residue.atom_count ||
            bond.atom2 < residue.atom_begin || bond.atom2 >= residue.atom_begin + residue.atom_count) {
            continue;
        }
        residue_type.add_connectivity(molecule.atoms[bond.atom1].name, molecule.atoms[bond.atom2].name);
    }
    return residue_type;
}

std::vector<std::string> split_ws(const std::string& line) {
    std::istringstream iss(line);
    std::vector<std::string> out;
    std::string word;
    while (iss >> word) {
        out.push_back(word);
    }
    return out;
}

std::string trim_copy(const std::string& input) {
    const auto first = input.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}

std::vector<std::string> split_dash_atoms(const std::string& input) {
    std::vector<std::string> atoms;
    std::string token;
    std::istringstream stream(input);
    while (std::getline(stream, token, '-')) {
        atoms.push_back(trim_copy(token));
    }
    return atoms;
}

std::pair<std::array<std::string, 4>, std::vector<std::string>> amber_atoms_words(
    const std::string& line, std::size_t atom_width, const std::array<std::string, 4>* last_atoms = nullptr) {
    std::array<std::string, 4> atoms{"", "", "", ""};
    std::vector<std::string> words;
    if (!line.empty() && line[0] == ' ') {
        if (last_atoms != nullptr) {
            atoms = *last_atoms;
        }
        words = split_ws(line.size() > atom_width ? line.substr(atom_width) : "");
        return {atoms, words};
    }
    const auto atom_words = split_dash_atoms(line.substr(0, std::min(atom_width, line.size())));
    for (std::size_t i = 0; i < std::min<std::size_t>(4, atom_words.size()); ++i) {
        atoms[i] = atom_words[i];
    }
    words = split_ws(line.size() > atom_width ? line.substr(atom_width) : "");
    return {atoms, words};
}

void upsert_bond_parameter(const std::string& atom_type1, const std::string& atom_type2, const BondTerm& term) {
    auto& parameters = bond_parameters();
    auto it = std::find_if(parameters.begin(), parameters.end(), [&](const BondParameter& parameter) {
        return (parameter.atom_type1 == atom_type1 && parameter.atom_type2 == atom_type2) ||
               (parameter.atom_type1 == atom_type2 && parameter.atom_type2 == atom_type1);
    });
    if (it == parameters.end()) {
        parameters.push_back({atom_type1, atom_type2, term, parameters.size() + 1});
    } else {
        it->term = term;
        it->order = parameters.size() + 1;
    }
}

void upsert_angle_parameter(const std::array<std::string, 3>& types, const AngleTerm& term) {
    auto& parameters = angle_parameters();
    auto reverse = std::array<std::string, 3>{types[2], types[1], types[0]};
    auto it = std::find_if(parameters.begin(), parameters.end(), [&](const AngleParameter& parameter) {
        return parameter.types == types || parameter.types == reverse;
    });
    if (it == parameters.end()) {
        parameters.push_back({types, term, parameters.size() + 1});
    } else {
        it->term = term;
        it->order = parameters.size() + 1;
    }
}

void upsert_proper_parameter(const std::array<std::string, 4>& types, const DihedralTerm& term, bool reset) {
    auto& parameters = proper_parameters();
    auto it = std::find_if(parameters.begin(), parameters.end(), [&](const ProperParameter& parameter) {
        return parameter.types == types;
    });
    if (it == parameters.end()) {
        parameters.push_back({types, {}, parameters.size() + 1});
        it = std::prev(parameters.end());
    }
    if (reset) {
        it->terms.clear();
    }
    it->terms.push_back(term);
}

void add_improper_parameter(const std::array<std::string, 4>& types, const DihedralTerm& term) {
    improper_parameters().push_back({types, term, improper_parameters().size() + 1});
}

void upsert_nb14_parameter(const std::string& atom_type1, const std::string& atom_type4, const NB14Scale& scale) {
    auto& parameters = nb14_parameters();
    auto it = std::find_if(parameters.begin(), parameters.end(), [&](const NB14Parameter& parameter) {
        return parameter.atom_type1 == atom_type1 && parameter.atom_type4 == atom_type4;
    });
    if (it == parameters.end()) {
        parameters.push_back({atom_type1, atom_type4, scale, parameters.size() + 1});
    } else {
        it->scale = scale;
    }
}

void upsert_lj_atom_type(const std::string& atom_type, const std::string& lj_type) {
    lj_type_by_atom_type()[atom_type] = lj_type;
}

void upsert_mass(const std::string& atom_type, double mass) {
    mass_by_atom_type()[atom_type] = mass;
}

void upsert_lj_parameter(const std::string& lj_type, double epsilon, double rmin) {
    lj_parameters()[lj_type] = {epsilon, rmin};
}

void upsert_amber_cmap_parameter(const std::string& residue, std::uint32_t resolution,
                                 const std::vector<double>& parameters) {
    amber_cmap_parameters().insert_or_assign(
        "C-N-" + residue + "@XC-C-N",
        AmberCMapParameter{resolution, parameters});
}

void upsert_amber_cmap_key(const std::string& key, std::uint32_t resolution,
                           const std::vector<double>& parameters) {
    amber_cmap_parameters().insert_or_assign(key, AmberCMapParameter{resolution, parameters});
}

void flush_amber_cmap_block(const std::vector<std::string>& residues, std::uint32_t resolution,
                            const std::vector<double>& parameters) {
    if (residues.empty() || parameters.empty()) {
        return;
    }
    for (const auto& residue : residues) {
        upsert_amber_cmap_parameter(residue, resolution, parameters);
    }
}

std::optional<NB14Scale> parse_nb14_scale(const std::string& line) {
    const auto scee_pos = line.find("SCEE=");
    const auto scnb_pos = line.find("SCNB=");
    if (scee_pos == std::string::npos && scnb_pos == std::string::npos) {
        return std::nullopt;
    }
    NB14Scale scale;
    if (scnb_pos != std::string::npos) {
        scale.k_lj = 1.0 / std::stod(line.substr(scnb_pos + 5));
    }
    if (scee_pos != std::string::npos) {
        scale.k_ee = 1.0 / std::stod(line.substr(scee_pos + 5));
    }
    return scale;
}

int wildcard_score(const std::array<std::string, 4>& parameter, const std::array<std::string, 4>& query) {
    int score = 0;
    for (std::size_t i = 0; i < 4; ++i) {
        if (parameter[i] == "X") {
            continue;
        }
        if (parameter[i] != query[i]) {
            return -1;
        }
        ++score;
    }
    return score;
}

int wildcard_pair_score(const std::string& parameter1, const std::string& parameter4, const std::string& query1,
                        const std::string& query4) {
    int score = 0;
    if (parameter1 != "X") {
        if (parameter1 != query1) {
            return -1;
        }
        ++score;
    }
    if (parameter4 != "X") {
        if (parameter4 != query4) {
            return -1;
        }
        ++score;
    }
    return score;
}

int wildcard_angle_score(const std::array<std::string, 3>& parameter, const std::array<std::string, 3>& query) {
    int score = 0;
    for (std::size_t i = 0; i < 3; ++i) {
        if (parameter[i] == "X") {
            continue;
        }
        if (parameter[i] != query[i]) {
            return -1;
        }
        ++score;
    }
    return score;
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

AtomId find_residue_atom_by_name(const Molecule& molecule, const Residue& residue, const std::string& name) {
    for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
        const AtomId atom_id = residue.atom_begin + local;
        if (molecule.atoms[atom_id].name == name) {
            return atom_id;
        }
    }
    return static_cast<AtomId>(molecule.atoms.size());
}

}  // namespace

void register_amber_parmdat_file(const std::filesystem::path& filename) {
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open Amber parmdat file: " + filename.string());
    }
    std::unique_lock lock(registry_mutex());

    std::string line;
    std::getline(input, line);
    while (std::getline(input, line) && !trim_copy(line).empty()) {
        const auto words = split_ws(line);
        if (!words.empty()) {
            upsert_lj_atom_type(words[0], words[0]);
        }
    }
    std::getline(input, line);
    while (std::getline(input, line) && !trim_copy(line).empty()) {
        auto [atoms, words] = amber_atoms_words(line, 5);
        if (words.size() >= 2) {
            upsert_bond_parameter(atoms[0], atoms[1], {std::stod(words[0]), std::stod(words[1])});
        }
    }
    while (std::getline(input, line) && !trim_copy(line).empty()) {
        auto [atoms4, words] = amber_atoms_words(line, 8);
        if (words.size() >= 2) {
            upsert_angle_parameter({atoms4[0], atoms4[1], atoms4[2]},
                                   {std::stod(words[0]), std::stod(words[1]) / 180.0 *
                                                            3.14159265358979323846});
        }
    }

    bool reset = true;
    std::array<std::string, 4> last_atoms{"", "", "", ""};
    while (std::getline(input, line)) {
        if (trim_copy(line).empty()) {
            break;
        }
        auto [atoms, words] = amber_atoms_words(line, 11, &last_atoms);
        last_atoms = atoms;
        if (words.size() < 4) {
            continue;
        }
        const auto divisor = std::stoi(words[0]);
        const auto k = std::stod(words[1]) / divisor;
        const auto phase = std::stod(words[2]) / 180.0 * 3.14159265358979323846;
        const auto periodicity = std::abs(static_cast<int>(std::stod(words[3])));
        if (const auto scale = parse_nb14_scale(line)) {
            upsert_nb14_parameter(atoms[0], atoms[3], *scale);
        }
        upsert_proper_parameter(atoms, {periodicity, k, phase}, reset);
        reset = std::stod(words[3]) >= 0.0;
    }

    while (std::getline(input, line)) {
        if (trim_copy(line).empty()) {
            break;
        }
        auto [atoms, words] = amber_atoms_words(line, 11);
        if (words.size() < 3) {
            continue;
        }
        add_improper_parameter(
            atoms, {static_cast<int>(std::stod(words[2])), std::stod(words[0]), std::stod(words[1]) / 180.0 *
                                                                       3.14159265358979323846});
    }

    std::getline(input, line);
    std::getline(input, line);
    while (std::getline(input, line) && !trim_copy(line).empty()) {
        auto words = split_ws(line);
        if (words.empty()) {
            continue;
        }
        const auto lj_type = words.front();
        for (std::size_t i = 1; i < words.size(); ++i) {
            upsert_lj_atom_type(words[i], lj_type);
        }
    }

    if (std::getline(input, line)) {
        const auto words = split_ws(line);
        if (words.size() >= 2 && words[1] != "RE") {
            throw std::runtime_error("unsupported Amber LJ parameter format in parmdat: " + words[1]);
        }
    }
    while (std::getline(input, line) && !trim_copy(line).empty()) {
        const auto words = split_ws(line);
        if (words.size() >= 3) {
            upsert_lj_parameter(words[0], std::stod(words[2]), std::stod(words[1]));
        }
    }
}

void register_amber_frcmod_file(const std::filesystem::path& filename) {
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open Amber frcmod file: " + filename.string());
    }
    std::unique_lock lock(registry_mutex());

    std::string line;
    std::string flag;
    std::string cmap_flag;
    std::vector<std::string> cmap_residues;
    std::uint32_t cmap_resolution = 24;
    std::vector<double> cmap_parameters;
    bool reset = true;
    std::getline(input, line);
    while (std::getline(input, line)) {
        const auto trimmed = trim_copy(line);
        if (trimmed.empty()) {
            continue;
        }
        const auto words0 = split_ws(line);
        if (flag.rfind("CMAP", 0) != 0 && words0.size() == 1) {
            flag = trimmed;
            continue;
        }
        if (flag.rfind("MASS", 0) == 0) {
            if (!words0.empty()) {
                upsert_lj_atom_type(words0[0], words0[0]);
                if (words0.size() >= 2) {
                    upsert_mass(words0[0], std::stod(words0[1]));
                }
            }
        } else if (flag.rfind("BOND", 0) == 0) {
            auto [atoms, words] = amber_atoms_words(line, 5);
            if (words.size() >= 2) {
                upsert_bond_parameter(atoms[0], atoms[1], {std::stod(words[0]), std::stod(words[1])});
            }
        } else if (flag.rfind("ANGL", 0) == 0) {
            auto [atoms4, words] = amber_atoms_words(line, 8);
            if (words.size() >= 2) {
                upsert_angle_parameter({atoms4[0], atoms4[1], atoms4[2]},
                                       {std::stod(words[0]), std::stod(words[1]) / 180.0 *
                                                                3.14159265358979323846});
            }
        } else if (flag.rfind("DIHE", 0) == 0) {
            auto [atoms, words] = amber_atoms_words(line, 11);
            if (words.size() < 4) {
                continue;
            }
            const auto divisor = std::stoi(words[0]);
            const auto k = std::stod(words[1]) / divisor;
            const auto phase = std::stod(words[2]) / 180.0 * 3.14159265358979323846;
            const auto periodicity = std::abs(static_cast<int>(std::stod(words[3])));
            if (const auto scale = parse_nb14_scale(line)) {
                upsert_nb14_parameter(atoms[0], atoms[3], *scale);
            }
            upsert_proper_parameter(atoms, {periodicity, k, phase}, reset);
            reset = std::stod(words[3]) >= 0.0;
        } else if (flag.rfind("IMPR", 0) == 0) {
            auto [atoms, words] = amber_atoms_words(line, 11);
            if (words.size() < 3) {
                continue;
            }
            add_improper_parameter(
                atoms, {static_cast<int>(std::stod(words[2])), std::stod(words[0]), std::stod(words[1]) / 180.0 *
                                                                           3.14159265358979323846});
        } else if (flag.rfind("NONB", 0) == 0) {
            if (words0.size() >= 3) {
                upsert_lj_atom_type(words0[0], words0[0]);
                upsert_lj_parameter(words0[0], std::stod(words0[2]), std::stod(words0[1]));
            }
        } else if (flag.rfind("CMAP", 0) == 0) {
            if (trimmed.rfind("%FLAG", 0) == 0) {
                if (trimmed.find("CMAP_COUNT") != std::string::npos) {
                    flush_amber_cmap_block(cmap_residues, cmap_resolution, cmap_parameters);
                    cmap_residues.clear();
                    cmap_parameters.clear();
                    cmap_resolution = 24;
                    cmap_flag = "CMAP_COUNT";
                } else if (trimmed.find("CMAP_RESLIST") != std::string::npos) {
                    cmap_flag = "CMAP_RESLIST";
                } else if (trimmed.find("CMAP_RESOLUTION") != std::string::npos) {
                    if (!words0.empty()) {
                        cmap_resolution = static_cast<std::uint32_t>(std::stoul(words0.back()));
                    }
                    cmap_flag = "CMAP_RESOLUTION";
                } else if (trimmed.find("CMAP_PARAMETER") != std::string::npos) {
                    cmap_flag = "CMAP_PARAMETER";
                } else if (trimmed.find("CMAP_TITLE") != std::string::npos) {
                    cmap_flag = "CMAP_TITLE";
                }
            } else if (cmap_flag == "CMAP_RESLIST") {
                cmap_residues.insert(cmap_residues.end(), words0.begin(), words0.end());
            } else if (cmap_flag == "CMAP_PARAMETER") {
                for (const auto& word : words0) {
                    cmap_parameters.push_back(std::stod(word));
                }
            }
        }
    }
    if (flag.rfind("CMAP", 0) == 0) {
        flush_amber_cmap_block(cmap_residues, cmap_resolution, cmap_parameters);
    }
}

void register_amber_lj_parameter(const std::string& atom_type, const std::string& lj_type, double epsilon,
                                 double rmin) {
    std::unique_lock lock(registry_mutex());
    upsert_lj_atom_type(atom_type, lj_type);
    upsert_lj_parameter(lj_type, epsilon, rmin);
}

void register_amber_bond_parameter(const std::string& atom_type1, const std::string& atom_type2, double k,
                                   double length) {
    std::unique_lock lock(registry_mutex());
    upsert_bond_parameter(atom_type1, atom_type2, {k, length});
}

void register_amber_cmap_parameter(const std::string& key, std::uint32_t resolution,
                                   const std::vector<double>& parameters) {
    if (resolution == 0 || parameters.size() != static_cast<std::size_t>(resolution) * resolution) {
        throw std::invalid_argument("Amber CMAP parameter count should equal resolution * resolution");
    }
    std::unique_lock lock(registry_mutex());
    upsert_amber_cmap_key(key, resolution, parameters);
}

bool has_amber_cmap_parameters() {
    std::shared_lock lock(registry_mutex());
    return !amber_cmap_parameters().empty();
}

void apply_amber_cmaps(Molecule& molecule) {
    if (!molecule.cmaps.empty()) {
        return;
    }
    std::shared_lock lock(registry_mutex());
    if (amber_cmap_parameters().empty() || molecule.residues.size() < 3) {
        return;
    }

    std::unordered_map<std::string, std::uint32_t> type_by_key;
    for (std::size_t residue_id = 1; residue_id + 1 < molecule.residues.size(); ++residue_id) {
        const auto& previous = molecule.residues[residue_id - 1];
        const auto& current = molecule.residues[residue_id];
        const auto& next = molecule.residues[residue_id + 1];
        const AtomId atom0 = find_residue_atom_by_name(molecule, previous, "C");
        const AtomId atom1 = find_residue_atom_by_name(molecule, current, "N");
        const AtomId atom2 = find_residue_atom_by_name(molecule, current, "CA");
        const AtomId atom3 = find_residue_atom_by_name(molecule, current, "C");
        const AtomId atom4 = find_residue_atom_by_name(molecule, next, "N");
        if (atom0 >= molecule.atoms.size() || atom1 >= molecule.atoms.size() || atom2 >= molecule.atoms.size() ||
            atom3 >= molecule.atoms.size() || atom4 >= molecule.atoms.size()) {
            continue;
        }
        const std::string key = molecule.atoms[atom0].type + "-" + molecule.atoms[atom1].type + "-" +
                                current.name + "@" + molecule.atoms[atom2].type + "-" +
                                molecule.atoms[atom3].type + "-" + molecule.atoms[atom4].type;
        const auto cmap_it = amber_cmap_parameters().find(key);
        if (cmap_it == amber_cmap_parameters().end()) {
            continue;
        }
        auto type_it = type_by_key.find(key);
        if (type_it == type_by_key.end()) {
            type_it = type_by_key.emplace(
                key, molecule.add_cmap_type(cmap_it->second.resolution, cmap_it->second.parameters)).first;
        }
        molecule.add_cmap(atom0, atom1, atom2, atom3, atom4, type_it->second);
    }
}

std::vector<DihedralTerm> find_amber_proper_terms(const std::array<std::string, 4>& atom_types) {
    std::shared_lock lock(registry_mutex());
    const std::array<std::string, 4> reverse{atom_types[3], atom_types[2], atom_types[1], atom_types[0]};
    int best_score = -1;
    std::size_t best_order = 0;
    std::vector<DihedralTerm> out;
    for (const auto& parameter : proper_parameters()) {
        const int score = std::max(wildcard_score(parameter.types, atom_types), wildcard_score(parameter.types, reverse));
        if (score < 0 || parameter.terms.empty()) {
            continue;
        }
        if (score > best_score || (score == best_score && parameter.order >= best_order)) {
            if (score > best_score || parameter.order > best_order) {
                out.clear();
            }
            best_score = score;
            best_order = parameter.order;
            out.insert(out.end(), parameter.terms.begin(), parameter.terms.end());
        }
    }
    return out;
}

std::optional<DihedralTerm> find_amber_improper_term(const std::array<std::string, 4>& atom_types) {
    std::shared_lock lock(registry_mutex());
    int best_score = -1;
    std::size_t best_order = 0;
    std::optional<DihedralTerm> out;
    for (const auto& parameter : improper_parameters()) {
        const int score = wildcard_score(parameter.types, atom_types);
        if (score < 0) {
            continue;
        }
        if (score > best_score || (score == best_score && parameter.order >= best_order)) {
            best_score = score;
            best_order = parameter.order;
            out = parameter.term;
        }
    }
    return out;
}

std::optional<AmberImproperMatch> find_amber_improper_match(const std::array<std::string, 4>& atom_types) {
    std::shared_lock lock(registry_mutex());
    for (auto it = improper_parameters().rbegin(); it != improper_parameters().rend(); ++it) {
        if (it->types == atom_types) {
            return AmberImproperMatch{it->term, true, 0};
        }
    }

    std::optional<AmberImproperMatch> out;
    for (std::uint32_t mask = 0; mask < 16; ++mask) {
        std::array<std::string, 4> candidate = atom_types;
        int x_count = 0;
        for (std::size_t i = 0; i < 4; ++i) {
            if ((mask & (1U << (3U - i))) != 0U) {
                candidate[i] = "X";
                ++x_count;
            }
        }
        for (auto it = improper_parameters().rbegin(); it != improper_parameters().rend(); ++it) {
            if (it->types == candidate) {
                return AmberImproperMatch{it->term, false, x_count};
            }
        }
    }
    return out;
}

std::optional<NB14Scale> find_amber_nb14_scale(const std::string& atom_type1, const std::string& atom_type4) {
    std::shared_lock lock(registry_mutex());
    int best_score = -1;
    std::size_t best_order = 0;
    std::optional<NB14Scale> out;
    for (const auto& parameter : nb14_parameters()) {
        const int score = std::max(wildcard_pair_score(parameter.atom_type1, parameter.atom_type4, atom_type1, atom_type4),
                                   wildcard_pair_score(parameter.atom_type1, parameter.atom_type4, atom_type4, atom_type1));
        if (score < 0) {
            continue;
        }
        if (score > best_score || (score == best_score && parameter.order >= best_order)) {
            best_score = score;
            best_order = parameter.order;
            out = parameter.scale;
        }
    }
    return out;
}

std::optional<BondTerm> find_amber_bond_term(const std::string& atom_type1, const std::string& atom_type2) {
    std::shared_lock lock(registry_mutex());
    int best_score = -1;
    std::size_t best_order = 0;
    std::optional<BondTerm> out;
    for (const auto& parameter : bond_parameters()) {
        const int score = std::max(wildcard_pair_score(parameter.atom_type1, parameter.atom_type2, atom_type1, atom_type2),
                                   wildcard_pair_score(parameter.atom_type1, parameter.atom_type2, atom_type2, atom_type1));
        if (score < 0) {
            continue;
        }
        if (score > best_score || (score == best_score && parameter.order >= best_order)) {
            best_score = score;
            best_order = parameter.order;
            out = parameter.term;
        }
    }
    return out;
}

std::optional<AngleTerm> find_amber_angle_term(const std::array<std::string, 3>& atom_types) {
    std::shared_lock lock(registry_mutex());
    const std::array<std::string, 3> reverse{atom_types[2], atom_types[1], atom_types[0]};
    int best_score = -1;
    std::size_t best_order = 0;
    std::optional<AngleTerm> out;
    for (const auto& parameter : angle_parameters()) {
        const int score = std::max(wildcard_angle_score(parameter.types, atom_types),
                                   wildcard_angle_score(parameter.types, reverse));
        if (score < 0) {
            continue;
        }
        if (score > best_score || (score == best_score && parameter.order >= best_order)) {
            best_score = score;
            best_order = parameter.order;
            out = parameter.term;
        }
    }
    return out;
}

std::string find_amber_lj_type(const std::string& atom_type) {
    std::shared_lock lock(registry_mutex());
    const auto it = lj_type_by_atom_type().find(atom_type);
    if (it != lj_type_by_atom_type().end()) {
        return it->second;
    }
    return atom_type;
}

std::optional<std::pair<double, double>> find_amber_lj_parameter(const std::string& lj_type) {
    std::shared_lock lock(registry_mutex());
    const auto it = lj_parameters().find(lj_type);
    if (it == lj_parameters().end()) {
        return std::nullopt;
    }
    return it->second;
}

std::optional<double> find_amber_atom_type_mass(const std::string& atom_type) {
    std::shared_lock lock(registry_mutex());
    const auto it = mass_by_atom_type().find(atom_type);
    if (it == mass_by_atom_type().end() || it->second <= 0.0) {
        return std::nullopt;
    }
    return it->second;
}

void register_ff14sb() {
    std::unique_lock lock(registry_mutex());
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
    std::unique_lock lock(registry_mutex());

    ResidueType wat("WAT");
    wat.add_atom("O", "OW", 0.0000, 0.0000, 0.0000, -0.834, 16.0);
    wat.add_atom("H1", "HW", 0.9572, 0.0000, 0.0000, 0.417, 1.008);
    wat.add_atom("H2", "HW", -0.239988, 0.926627, 0.0000, 0.417, 1.008);
    wat.add_connectivity("O", "H1");
    wat.add_connectivity("O", "H2");
    wat.add_connectivity("H1", "H2");
    put_template(std::move(wat));
    upsert_lj_atom_type("OW", "OW");
    upsert_lj_atom_type("HW", "HW");
    upsert_lj_parameter("OW", 0.152, 1.7683);
    upsert_lj_parameter("HW", 0.0, 0.0);

    ResidueType na("NA");
    na.add_atom("NA", "Na+", 0.0, 0.0, 0.0, 1.0, 22.99);
    put_template(std::move(na));
    upsert_lj_atom_type("Na+", "Na+");

    ResidueType cl("CL");
    cl.add_atom("CL", "Cl-", 0.0, 0.0, 0.0, -1.0, 35.45);
    put_template(std::move(cl));
    upsert_lj_atom_type("Cl-", "Cl-");
}

void register_residue_templates_from_mol2_text(const std::string& text) {
    const auto molecule = load_mol2_text(text);
    std::vector<std::pair<ResidueId, ResidueType>> residue_types;
    residue_types.reserve(molecule.residues.size());
    std::unordered_map<ResidueId, std::size_t> residue_to_type;

    for (ResidueId residue_id = 0; residue_id < molecule.residues.size(); ++residue_id) {
        const auto& residue = molecule.residues[residue_id];
        residue_to_type[residue_id] = residue_types.size();
        residue_types.emplace_back(residue_id, ResidueType(residue.name));
        auto& residue_type = residue_types.back().second;
        for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
            const auto& atom = molecule.atoms[residue.atom_begin + local];
            residue_type.add_atom(atom.name, atom.type, atom.x, atom.y, atom.z, atom.charge, atom.mass);
        }
    }

    for (const auto& bond : molecule.explicit_bonds) {
        const auto res1 = molecule.atoms[bond.atom1].residue;
        const auto res2 = molecule.atoms[bond.atom2].residue;
        if (res1 != res2) {
            continue;
        }
        auto& residue_type = residue_types[residue_to_type.at(res1)].second;
        try {
            residue_type.add_connectivity(molecule.atoms[bond.atom1].name, molecule.atoms[bond.atom2].name);
        } catch (const std::exception&) {
        }
    }

    std::unique_lock lock(registry_mutex());
    for (auto& residue_type : residue_types) {
        put_template(std::move(residue_type.second));
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

void register_template_molecule_from_mol2_file(const std::filesystem::path& filename) {
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open mol2 template file: " + filename.string());
    }
    std::ostringstream buffer;
    buffer << input.rdbuf();
    auto molecule = load_mol2_text(buffer.str());
    if (molecule.residues.empty()) {
        throw std::invalid_argument("mol2 template file contains no residues: " + filename.string());
    }

    std::vector<ResidueType> residue_types;
    residue_types.reserve(molecule.residues.size());
    for (const auto& residue : molecule.residues) {
        residue_types.push_back(residue_type_from_molecule_residue(molecule, residue));
    }

    std::unique_lock lock(registry_mutex());
    for (auto& residue_type : residue_types) {
        put_template(std::move(residue_type));
    }
    molecule_templates().insert_or_assign(molecule.residues.front().name, std::move(molecule));
}

void register_template_virtual_atom2(const std::string& template_name, const std::string& virtual_atom,
                                     const std::string& atom0, const std::string& atom1, const std::string& atom2,
                                     double k1, double k2) {
    std::unique_lock lock(registry_mutex());
    auto it = molecule_templates().find(template_name);
    if (it == molecule_templates().end()) {
        throw std::out_of_range("molecule template not found: " + template_name);
    }
    auto& molecule = it->second;
    if (molecule.residues.empty()) {
        throw std::invalid_argument("molecule template contains no residues: " + template_name);
    }
    const auto& residue = molecule.residues.front();
    const auto find_atom_by_name = [&](const std::string& name) {
        for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
            const AtomId atom_id = residue.atom_begin + local;
            if (molecule.atoms[atom_id].name == name) {
                return atom_id;
            }
        }
        throw std::out_of_range("atom name not found in molecule template " + template_name + ": " + name);
    };
    molecule.add_virtual_atom2(find_atom_by_name(virtual_atom), find_atom_by_name(atom0),
                               find_atom_by_name(atom1), find_atom_by_name(atom2), k1, k2);
}

void configure_residue_template_head(const std::string& template_name, const std::string& atom,
                                     double length, const std::string& next) {
    std::unique_lock lock(registry_mutex());
    auto it = templates().find(template_name);
    if (it == templates().end()) {
        throw std::out_of_range("residue template not found: " + template_name);
    }
    it->second.set_head(atom, length, next);
}

void configure_residue_template_tail(const std::string& template_name, const std::string& atom,
                                     double length, const std::string& next) {
    std::unique_lock lock(registry_mutex());
    auto it = templates().find(template_name);
    if (it == templates().end()) {
        throw std::out_of_range("residue template not found: " + template_name);
    }
    it->second.set_tail(atom, length, next);
}

void configure_residue_template_connect_atom(const std::string& template_name, const std::string& key,
                                             const std::string& atom) {
    std::unique_lock lock(registry_mutex());
    auto it = templates().find(template_name);
    if (it == templates().end()) {
        throw std::out_of_range("residue template not found: " + template_name);
    }
    it->second.set_connect_atom(key, atom);
}

void register_residue_template_alias(const std::string& alias_name, const std::string& template_name) {
    std::unique_lock lock(registry_mutex());
    const auto it = templates().find(template_name);
    if (it == templates().end()) {
        throw std::out_of_range("residue template not found: " + template_name);
    }
    ResidueType alias(alias_name);
    for (const auto& atom : it->second.atoms()) {
        alias.add_atom(atom.name, atom.type, atom.x, atom.y, atom.z, atom.charge, atom.mass);
    }
    for (const auto& bond : it->second.bonds()) {
        alias.add_connectivity(it->second.atoms()[bond.atom1].name, it->second.atoms()[bond.atom2].name);
    }
    if (!it->second.head().empty()) {
        alias.set_head(it->second.head(), it->second.head_length(), it->second.head_next());
    }
    if (!it->second.tail().empty()) {
        alias.set_tail(it->second.tail(), it->second.tail_length(), it->second.tail_next());
    }
    for (const auto& [key, atom] : it->second.connect_atoms()) {
        alias.set_connect_atom(key, atom);
    }
    templates().insert_or_assign(alias_name, std::move(alias));

    const auto molecule_it = molecule_templates().find(template_name);
    if (molecule_it != molecule_templates().end()) {
        Molecule alias_molecule = molecule_it->second;
        alias_molecule.name = alias_name;
        for (auto& residue : alias_molecule.residues) {
            residue.name = alias_name;
            residue.type_name = alias_name;
            residue.original_name = alias_name;
        }
        molecule_templates().insert_or_assign(alias_name, std::move(alias_molecule));
    }
}

bool has_template(const std::string& name) {
    std::shared_lock lock(registry_mutex());
    return templates().find(name) != templates().end();
}

std::size_t template_atom_count(const std::string& name) {
    std::shared_lock lock(registry_mutex());
    const auto it = templates().find(name);
    if (it == templates().end()) {
        throw std::out_of_range("residue template not found: " + name);
    }
    return it->second.atom_count();
}

std::vector<std::string> registered_template_names() {
    std::shared_lock lock(registry_mutex());
    std::vector<std::string> names;
    names.reserve(templates().size());
    for (const auto& [name, _template] : templates()) {
        names.push_back(name);
    }
    std::sort(names.begin(), names.end());
    return names;
}

const ResidueType& get_residue_template(const std::string& name) {
    std::shared_lock lock(registry_mutex());
    const auto it = templates().find(name);
    if (it == templates().end()) {
        throw std::out_of_range("residue template not found: " + name);
    }
    return it->second;
}

Molecule get_template_molecule(const std::string& name) {
    {
        std::shared_lock lock(registry_mutex());
        const auto molecule_it = molecule_templates().find(name);
        if (molecule_it != molecule_templates().end()) {
            return molecule_it->second;
        }
    }
    const auto& residue_type = get_residue_template(name);
    Molecule molecule(name);
    molecule.append_residue_from_type(residue_type, 0.0, 0.0, 0.0);
    if (!molecule.atoms.empty()) {
        molecule.set_box_padding(0.0, false);
    }
    return molecule;
}

}  // namespace xpongecpp
