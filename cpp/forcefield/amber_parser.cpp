#include "amber_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <mutex>
#include <numeric>
#include <optional>
#include <shared_mutex>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace xpongecpp {
namespace {

std::vector<std::string> split_ws(const std::string& line) {
    std::istringstream iss(line);
    std::vector<std::string> out;
    std::string word;
    while (iss >> word) {
        out.push_back(word);
    }
    return out;
}

LJCombiningRule& current_lj_combining_rule() {
    static LJCombiningRule rule = LJCombiningRule::LorentzBerthelot;
    return rule;
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
    const auto atom_field = line.substr(0, std::min(atom_width, line.size()));
    if (trim_copy(atom_field).empty()) {
        if (last_atoms == nullptr) {
            throw std::runtime_error("Amber parameter continuation line has no preceding atom types");
        }
        atoms = *last_atoms;
        words = split_ws(line.size() > atom_width ? line.substr(atom_width) : "");
        return {atoms, words};
    }
    const auto atom_words = split_dash_atoms(atom_field);
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
    const std::size_t next_order = std::accumulate(
        parameters.begin(), parameters.end(), std::size_t{0},
        [](std::size_t current, const ProperParameter& parameter) {
            return std::max(current, parameter.order);
        }) + 1;
    auto it = std::find_if(parameters.begin(), parameters.end(), [&](const ProperParameter& parameter) {
        return parameter.types == types;
    });
    if (it == parameters.end()) {
        parameters.push_back({types, {}, next_order});
        it = std::prev(parameters.end());
    }
    if (reset) {
        it->terms.clear();
        it->order = next_order;
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

void upsert_amber_cmap_parameter(const std::string& residue, std::uint32_t resolution,
                                 const std::vector<double>& parameters) {
    amber_cmap_parameters().insert_or_assign(
        "C-N-" + residue + "@XC-C-N",
        AmberCMapParameter{resolution, parameters});
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

void upsert_lj_atom_type(const std::string& atom_type, const std::string& lj_type) {
    lj_type_by_atom_type()[atom_type] = lj_type;
}

void upsert_mass(const std::string& atom_type, double mass) {
    mass_by_atom_type()[atom_type] = mass;
}

void upsert_lj_parameter(const std::string& lj_type, double epsilon, double rmin) {
    lj_parameters()[lj_type] = {epsilon, rmin};
}

void upsert_amber_cmap_key(const std::string& key, std::uint32_t resolution,
                           const std::vector<double>& parameters) {
    amber_cmap_parameters().insert_or_assign(key, AmberCMapParameter{resolution, parameters});
}

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
    std::optional<std::array<std::string, 4>> last_dihedral_atoms;
    std::getline(input, line);
    while (std::getline(input, line)) {
        const auto trimmed = trim_copy(line);
        if (trimmed.empty()) {
            continue;
        }
        const auto words0 = split_ws(line);
        if (flag.rfind("CMAP", 0) != 0 && words0.size() == 1) {
            flag = trimmed;
            if (flag.rfind("DIHE", 0) == 0) {
                last_dihedral_atoms.reset();
                reset = true;
            }
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
            auto [atoms, words] = amber_atoms_words(
                line, 11, last_dihedral_atoms ? &*last_dihedral_atoms : nullptr);
            last_dihedral_atoms = atoms;
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

void register_amber_angle_parameter(const std::array<std::string, 3>& atom_types, double k, double theta) {
    std::unique_lock lock(registry_mutex());
    upsert_angle_parameter(atom_types, {k, theta});
}

void register_amber_proper_dihedral_parameter(const std::array<std::string, 4>& atom_types, int periodicity,
                                              double k, double phase, bool reset) {
    std::unique_lock lock(registry_mutex());
    upsert_proper_parameter(atom_types, {periodicity, k, phase}, reset);
}

void register_amber_improper_dihedral_parameter(const std::array<std::string, 4>& atom_types, int periodicity,
                                                double k, double phase) {
    std::unique_lock lock(registry_mutex());
    add_improper_parameter(atom_types, {periodicity, k, phase});
}

void register_amber_nb14_scale(const std::string& atom_type1, const std::string& atom_type4, double k_lj,
                               double k_ee) {
    std::unique_lock lock(registry_mutex());
    upsert_nb14_parameter(atom_type1, atom_type4, {k_lj, k_ee});
}

void register_amber_cmap_parameter(const std::string& key, std::uint32_t resolution,
                                   const std::vector<double>& parameters) {
    if (resolution == 0 || parameters.size() != static_cast<std::size_t>(resolution) * resolution) {
        throw std::invalid_argument("Amber CMAP parameter count should equal resolution * resolution");
    }
    std::unique_lock lock(registry_mutex());
    upsert_amber_cmap_key(key, resolution, parameters);
}

void clear_amber_dihedral_parameters() {
    std::unique_lock lock(registry_mutex());
    proper_parameters().clear();
}

void clear_amber_improper_parameters() {
    std::unique_lock lock(registry_mutex());
    improper_parameters().clear();
}

void set_lj_combining_rule(LJCombiningRule rule) {
    current_lj_combining_rule() = rule;
}

LJCombiningRule lj_combining_rule() {
    return current_lj_combining_rule();
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
    return std::nullopt;
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

}  // namespace xpongecpp
