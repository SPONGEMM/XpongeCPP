#include "core.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
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

std::unordered_map<std::string, std::string>& lj_type_by_atom_type() {
    static std::unordered_map<std::string, std::string> parameters;
    return parameters;
}

std::unordered_map<std::string, std::pair<double, double>>& lj_parameters() {
    static std::unordered_map<std::string, std::pair<double, double>> parameters;
    return parameters;
}

void put_template(ResidueType residue_type) {
    templates().insert_or_assign(residue_type.name(), std::move(residue_type));
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

std::string proper_key(const std::array<std::string, 4>& types) {
    return types[0] + "-" + types[1] + "-" + types[2] + "-" + types[3];
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

void upsert_lj_parameter(const std::string& lj_type, double epsilon, double rmin) {
    lj_parameters()[lj_type] = {epsilon, rmin};
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

}  // namespace

void register_amber_parmdat_file(const std::filesystem::path& filename) {
    std::ifstream input(filename);
    if (!input) {
        throw std::runtime_error("failed to open Amber parmdat file: " + filename.string());
    }
    std::scoped_lock lock(registry_mutex());

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
    std::scoped_lock lock(registry_mutex());

    std::string line;
    std::string flag;
    bool reset = true;
    std::getline(input, line);
    while (std::getline(input, line)) {
        const auto trimmed = trim_copy(line);
        if (trimmed.empty()) {
            continue;
        }
        const auto words0 = split_ws(line);
        if (words0.size() == 1) {
            flag = trimmed;
            continue;
        }
        if (flag.rfind("MASS", 0) == 0) {
            if (!words0.empty()) {
                upsert_lj_atom_type(words0[0], words0[0]);
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
        }
    }
}

std::vector<DihedralTerm> find_amber_proper_terms(const std::array<std::string, 4>& atom_types) {
    std::scoped_lock lock(registry_mutex());
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
    std::scoped_lock lock(registry_mutex());
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
    std::scoped_lock lock(registry_mutex());
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
    std::scoped_lock lock(registry_mutex());
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
    std::scoped_lock lock(registry_mutex());
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
    std::scoped_lock lock(registry_mutex());
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
    std::scoped_lock lock(registry_mutex());
    const auto it = lj_type_by_atom_type().find(atom_type);
    if (it != lj_type_by_atom_type().end()) {
        return it->second;
    }
    return atom_type;
}

std::optional<std::pair<double, double>> find_amber_lj_parameter(const std::string& lj_type) {
    std::scoped_lock lock(registry_mutex());
    const auto it = lj_parameters().find(lj_type);
    if (it == lj_parameters().end()) {
        return std::nullopt;
    }
    return it->second;
}

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
    struct Mol2AtomInfo {
        std::string name;
        std::string type;
        std::string residue_name;
        std::string residue_key;
        double x{0.0};
        double y{0.0};
        double z{0.0};
        double charge{0.0};
    };
    std::istringstream input(text);
    std::string line;
    std::string section;
    std::unordered_map<int, Mol2AtomInfo> atoms_by_id;
    std::vector<std::array<int, 2>> bonds_by_id;

    while (std::getline(input, line)) {
        if (line.rfind("@<TRIPOS>", 0) == 0) {
            section = trim_copy(line.substr(9));
            continue;
        }
        const auto words = split_ws(line);
        if (words.empty()) {
            continue;
        }
        if (section == "ATOM") {
            if (words.size() < 9) {
                continue;
            }
            Mol2AtomInfo atom;
            atom.name = words[1];
            atom.x = std::stod(words[2]);
            atom.y = std::stod(words[3]);
            atom.z = std::stod(words[4]);
            atom.type = words[5];
            atom.residue_name = words[7];
            atom.residue_key = words[6] + ":" + words[7];
            atom.charge = std::stod(words[8]);
            atoms_by_id[std::stoi(words[0])] = std::move(atom);
        } else if (section == "BOND") {
            if (words.size() < 4) {
                continue;
            }
            bonds_by_id.push_back({std::stoi(words[1]), std::stoi(words[2])});
        }
    }

    std::vector<int> atom_ids;
    atom_ids.reserve(atoms_by_id.size());
    for (const auto& [id, atom] : atoms_by_id) {
        (void)atom;
        atom_ids.push_back(id);
    }
    std::sort(atom_ids.begin(), atom_ids.end());
    std::vector<std::string> residue_order;
    std::unordered_map<std::string, ResidueType> residue_types;
    std::unordered_map<int, std::string> atom_residue;
    std::unordered_map<int, std::string> atom_name;
    for (const auto id : atom_ids) {
        const auto& atom = atoms_by_id.at(id);
        if (residue_types.find(atom.residue_key) == residue_types.end()) {
            residue_order.push_back(atom.residue_key);
            residue_types.emplace(atom.residue_key, ResidueType(atom.residue_name));
        }
        const bool atomic_ion = atom.residue_name == "NA" || atom.residue_name == "CL";
        residue_types.at(atom.residue_key)
            .add_atom(atom.name, atom.type, atom.x, atom.y, atom.z, atom.charge,
                      default_mass_for_element(guess_element(atom.name, atomic_ion ? atom.type : "")));
        atom_residue[id] = atom.residue_key;
        atom_name[id] = atom.name;
    }
    for (const auto& bond : bonds_by_id) {
        const auto res1 = atom_residue.find(bond[0]);
        const auto res2 = atom_residue.find(bond[1]);
        if (res1 == atom_residue.end() || res2 == atom_residue.end() || res1->second != res2->second) {
            continue;
        }
        try {
            residue_types.at(res1->second).add_connectivity(atom_name.at(bond[0]), atom_name.at(bond[1]));
        } catch (const std::exception&) {
        }
    }

    std::scoped_lock lock(registry_mutex());
    for (const auto& residue_name : residue_order) {
        put_template(std::move(residue_types.at(residue_name)));
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
