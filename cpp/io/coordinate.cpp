#include "core.hpp"

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace xpongecpp {
namespace {

std::vector<double> parse_doubles(const std::string& line) {
    std::istringstream input(line);
    std::vector<double> values;
    double value = 0.0;
    while (input >> value) {
        values.push_back(value);
    }
    return values;
}

void apply_coordinate_data(const CoordinateData& data, Molecule& molecule) {
    if (data.coordinates.size() != molecule.atoms.size()) {
        throw std::invalid_argument("coordinate atom count does not match molecule atom count");
    }
    for (std::size_t i = 0; i < data.coordinates.size(); ++i) {
        molecule.atoms[i].x = data.coordinates[i][0];
        molecule.atoms[i].y = data.coordinates[i][1];
        molecule.atoms[i].z = data.coordinates[i][2];
    }
    if (data.has_box) {
        molecule.box_length = {data.box[0], data.box[1], data.box[2]};
        molecule.box_angle = {data.box[3], data.box[4], data.box[5]};
        molecule.has_box = true;
    }
}

}  // namespace

CoordinateData load_coordinate_text(const std::string& text) {
    std::istringstream input(text);
    std::string line;
    if (!std::getline(input, line)) {
        throw std::invalid_argument("empty coordinate file");
    }
    std::istringstream header(line);
    std::size_t atom_count = 0;
    if (!(header >> atom_count)) {
        throw std::invalid_argument("invalid coordinate atom count");
    }

    CoordinateData data;
    data.coordinates.reserve(atom_count);
    for (std::size_t i = 0; i < atom_count; ++i) {
        if (!std::getline(input, line)) {
            throw std::invalid_argument("coordinate file ended before all atoms were read");
        }
        const auto values = parse_doubles(line);
        if (values.size() < 3) {
            throw std::invalid_argument("coordinate line should contain at least 3 values");
        }
        data.coordinates.push_back({values[0], values[1], values[2]});
    }

    if (!std::getline(input, line)) {
        throw std::invalid_argument("coordinate file is missing box line");
    }
    const auto box = parse_doubles(line);
    if (box.size() != 3 && box.size() != 6) {
        throw std::invalid_argument("coordinate box line should contain 3 or 6 values");
    }
    data.box = {box[0], box[1], box[2], 90.0, 90.0, 90.0};
    data.has_box = true;
    if (box.size() == 6) {
        data.box = {box[0], box[1], box[2], box[3], box[4], box[5]};
    }
    return data;
}

CoordinateData load_coordinate_text(const std::string& text, Molecule& molecule) {
    auto data = load_coordinate_text(text);
    apply_coordinate_data(data, molecule);
    return data;
}

CoordinateData load_rst7_text(const std::string& text) {
    std::istringstream input(text);
    std::string line;
    if (!std::getline(input, line) || !std::getline(input, line)) {
        throw std::invalid_argument("rst7 file is missing header");
    }
    std::istringstream header(line);
    std::size_t atom_count = 0;
    if (!(header >> atom_count)) {
        throw std::invalid_argument("invalid rst7 atom count");
    }

    std::vector<double> values;
    values.reserve(atom_count * 3 + 6);
    while (std::getline(input, line)) {
        const auto parsed = parse_doubles(line);
        values.insert(values.end(), parsed.begin(), parsed.end());
    }
    if (values.size() < atom_count * 3) {
        throw std::invalid_argument("rst7 file ended before all coordinates were read");
    }

    CoordinateData data;
    data.coordinates.reserve(atom_count);
    for (std::size_t i = 0; i < atom_count; ++i) {
        data.coordinates.push_back({values[i * 3], values[i * 3 + 1], values[i * 3 + 2]});
    }
    const std::size_t box_begin = atom_count * 3;
    if (values.size() >= box_begin + 6) {
        data.box = {values[box_begin], values[box_begin + 1], values[box_begin + 2],
                    values[box_begin + 3], values[box_begin + 4], values[box_begin + 5]};
        data.has_box = true;
    } else if (values.size() >= box_begin + 3) {
        data.box = {values[box_begin], values[box_begin + 1], values[box_begin + 2], 90.0, 90.0, 90.0};
        data.has_box = true;
    }
    return data;
}

CoordinateData load_rst7_text(const std::string& text, Molecule& molecule) {
    auto data = load_rst7_text(text);
    apply_coordinate_data(data, molecule);
    return data;
}

}  // namespace xpongecpp
