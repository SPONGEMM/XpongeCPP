#include "core.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace xpongecpp {
namespace {

std::vector<std::string> split_ws(const std::string& line) {
    std::istringstream input(line);
    std::vector<std::string> words;
    std::string word;
    while (input >> word) {
        words.push_back(word);
    }
    return words;
}

double dot(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double norm(const std::array<double, 3>& a) {
    return std::sqrt(dot(a, a));
}

double angle_degrees(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    const double denom = norm(a) * norm(b);
    if (denom == 0.0) {
        return 90.0;
    }
    const double cosine = std::clamp(dot(a, b) / denom, -1.0, 1.0);
    return std::acos(cosine) * 180.0 / std::acos(-1.0);
}

std::array<double, 6> box_from_gro_words(const std::vector<std::string>& words, bool read_box_angle) {
    if (words.size() == 3) {
        return {std::stod(words[0]) * 10.0, std::stod(words[1]) * 10.0, std::stod(words[2]) * 10.0,
                90.0, 90.0, 90.0};
    }
    if (words.size() != 9) {
        throw std::runtime_error("unsupported GRO box format");
    }
    std::array<double, 9> vals{};
    for (std::size_t i = 0; i < vals.size(); ++i) {
        vals[i] = std::stod(words[i]) * 10.0;
    }
    const std::array<double, 3> v1{vals[0], vals[3], vals[4]};
    const std::array<double, 3> v2{vals[5], vals[1], vals[6]};
    const std::array<double, 3> v3{vals[7], vals[8], vals[2]};
    std::array<double, 6> box{norm(v1), norm(v2), norm(v3), 90.0, 90.0, 90.0};
    if (read_box_angle) {
        box[3] = angle_degrees(v2, v3);
        box[4] = angle_degrees(v1, v3);
        box[5] = angle_degrees(v1, v2);
    }
    return box;
}

std::array<double, 9> gro_words_from_box(const std::array<double, 3>& length,
                                         const std::array<double, 3>& angle) {
    const double pi = std::acos(-1.0);
    const double alpha = angle[0] * pi / 180.0;
    const double beta = angle[1] * pi / 180.0;
    const double gamma = angle[2] * pi / 180.0;
    const double a = length[0];
    const double b = length[1];
    const double c = length[2];
    const double sin_gamma = std::sin(gamma);
    if (std::abs(sin_gamma) < 1e-12) {
        throw std::invalid_argument("invalid box gamma angle for GRO output");
    }
    const std::array<double, 3> v1{a, 0.0, 0.0};
    const std::array<double, 3> v2{b * std::cos(gamma), b * sin_gamma, 0.0};
    const double v3x = c * std::cos(beta);
    const double v3y = c * (std::cos(alpha) - std::cos(beta) * std::cos(gamma)) / sin_gamma;
    const double v3z2 = c * c - v3x * v3x - v3y * v3y;
    const std::array<double, 3> v3{v3x, v3y, std::sqrt(std::max(0.0, v3z2))};
    return {v1[0] / 10.0, v2[1] / 10.0, v3[2] / 10.0,
            v1[1] / 10.0, v1[2] / 10.0, v2[0] / 10.0,
            v2[2] / 10.0, v3[0] / 10.0, v3[1] / 10.0};
}

double gro_coordinate_field(const std::string& line, std::size_t begin) {
    if (line.size() >= begin + 8) {
        return std::stod(line.substr(begin, 8)) * 10.0;
    }
    const auto words = split_ws(line);
    if (words.size() < 3) {
        throw std::runtime_error("malformed GRO atom line");
    }
    return std::stod(words[words.size() - 3 + (begin - 20) / 8]) * 10.0;
}

std::array<double, 3> molecule_min(const Molecule& molecule) {
    std::array<double, 3> minv{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
    };
    for (const auto& atom : molecule.atoms) {
        minv[0] = std::min(minv[0], atom.x);
        minv[1] = std::min(minv[1], atom.y);
        minv[2] = std::min(minv[2], atom.z);
    }
    return minv;
}

std::array<double, 3> molecule_max(const Molecule& molecule) {
    std::array<double, 3> maxv{
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity(),
    };
    for (const auto& atom : molecule.atoms) {
        maxv[0] = std::max(maxv[0], atom.x);
        maxv[1] = std::max(maxv[1], atom.y);
        maxv[2] = std::max(maxv[2], atom.z);
    }
    return maxv;
}

}  // namespace

GroData load_gro_text(const std::string& text) {
    return load_gro_text(text, true);
}

GroData load_gro_text(const std::string& text, bool read_box_angle) {
    std::istringstream input(text);
    std::string line;
    if (!std::getline(input, line)) {
        throw std::runtime_error("empty GRO input");
    }
    if (!std::getline(input, line)) {
        throw std::runtime_error("missing GRO atom count");
    }
    const std::size_t atom_count = static_cast<std::size_t>(std::stoul(split_ws(line).at(0)));
    GroData data;
    data.coordinates.reserve(atom_count);
    for (std::size_t i = 0; i < atom_count; ++i) {
        if (!std::getline(input, line)) {
            throw std::runtime_error("unexpected EOF in GRO atom block");
        }
        data.coordinates.push_back({gro_coordinate_field(line, 20), gro_coordinate_field(line, 28),
                                    gro_coordinate_field(line, 36)});
    }
    if (!std::getline(input, line)) {
        throw std::runtime_error("missing GRO box line");
    }
    data.box = box_from_gro_words(split_ws(line), read_box_angle);
    return data;
}

GroData load_gro_text(const std::string& text, Molecule& molecule, bool read_box_angle) {
    auto data = load_gro_text(text, read_box_angle);
    if (data.coordinates.size() != molecule.atoms.size()) {
        throw std::runtime_error("GRO atom count does not match molecule atom count");
    }
    for (std::size_t i = 0; i < molecule.atoms.size(); ++i) {
        molecule.atoms[i].x = data.coordinates[i][0];
        molecule.atoms[i].y = data.coordinates[i][1];
        molecule.atoms[i].z = data.coordinates[i][2];
    }
    molecule.box_length = {data.box[0], data.box[1], data.box[2]};
    molecule.box_angle = {data.box[3], data.box[4], data.box[5]};
    molecule.has_box = true;
    return data;
}

void save_gro(const Molecule& molecule, const std::filesystem::path& filename) {
    if (!molecule.validate()) {
        throw std::invalid_argument("cannot export invalid molecule");
    }
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("failed to open GRO output: " + filename.string());
    }
    out << "Generated By Xponge\n";
    out << molecule.atoms.size() << "\n";
    out << std::fixed << std::setprecision(3);
    std::array<double, 3> coordinate_shift{0.0, 0.0, 0.0};
    std::array<double, 3> output_box = molecule.box_length;
    const bool has_explicit_box = molecule.has_box || molecule.box_length[0] != 0.0 ||
                                  molecule.box_length[1] != 0.0 || molecule.box_length[2] != 0.0;
    if (!has_explicit_box) {
        const auto minv = molecule_min(molecule);
        const auto maxv = molecule_max(molecule);
        constexpr double boxspace = 3.0;
        coordinate_shift = {boxspace - minv[0], boxspace - minv[1], boxspace - minv[2]};
        output_box = {maxv[0] - minv[0] + boxspace * 2.0,
                      maxv[1] - minv[1] + boxspace * 2.0,
                      maxv[2] - minv[2] + boxspace * 2.0};
    }
    for (std::size_t residue_index = 0; residue_index < molecule.residues.size(); ++residue_index) {
        const auto& residue = molecule.residues[residue_index];
        for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
            const AtomId atom_id = residue.atom_begin + local;
            const auto& atom = molecule.atoms[atom_id];
            out << std::setw(5) << ((residue_index + 1) % 100000)
                << std::left << std::setw(5) << residue.name.substr(0, 5)
                << std::right << std::setw(5) << atom.name.substr(0, 5)
                << std::setw(5) << ((atom_id + 1) % 100000)
                << std::setw(8) << (atom.x + coordinate_shift[0]) / 10.0
                << std::setw(8) << (atom.y + coordinate_shift[1]) / 10.0
                << std::setw(8) << (atom.z + coordinate_shift[2]) / 10.0 << "\n";
        }
    }
    out << std::setprecision(5);
    const bool orthogonal = std::abs(molecule.box_angle[0] - 90.0) < 1e-8 &&
                            std::abs(molecule.box_angle[1] - 90.0) < 1e-8 &&
                            std::abs(molecule.box_angle[2] - 90.0) < 1e-8;
    if (orthogonal) {
        out << std::setw(10) << output_box[0] / 10.0
            << std::setw(10) << output_box[1] / 10.0
            << std::setw(10) << output_box[2] / 10.0 << "\n";
    } else {
        const auto vals = gro_words_from_box(output_box, molecule.box_angle);
        for (const double value : vals) {
            out << std::setw(10) << value;
        }
        out << "\n";
    }
}

}  // namespace xpongecpp
