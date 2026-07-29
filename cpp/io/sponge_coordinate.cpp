#include "sponge_writers.hpp"

#include <algorithm>
#include <limits>

namespace xpongecpp {

std::array<double, 6> sponge_coordinate_box_for_export(const Molecule& molecule) {
    std::array<double, 6> box{molecule.box_length[0], molecule.box_length[1], molecule.box_length[2],
                              molecule.box_angle[0], molecule.box_angle[1], molecule.box_angle[2]};
    if (molecule.has_box || molecule.atoms.empty()) {
        return box;
    }

    std::array<double, 3> minv{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::infinity()};
    std::array<double, 3> maxv{-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
                               -std::numeric_limits<double>::infinity()};
    for (const auto& atom : molecule.atoms) {
        minv[0] = std::min(minv[0], atom.x);
        minv[1] = std::min(minv[1], atom.y);
        minv[2] = std::min(minv[2], atom.z);
        maxv[0] = std::max(maxv[0], atom.x);
        maxv[1] = std::max(maxv[1], atom.y);
        maxv[2] = std::max(maxv[2], atom.z);
    }
    box[0] = maxv[0] - minv[0] + 6.0;
    box[1] = maxv[1] - minv[1] + 6.0;
    box[2] = maxv[2] - minv[2] + 6.0;
    return box;
}

std::array<double, 3> sponge_coordinate_shift_for_export(const Molecule& molecule) {
    if (molecule.has_box_origin) {
        return {-molecule.box_origin[0], -molecule.box_origin[1], -molecule.box_origin[2]};
    }
    if (molecule.has_box || molecule.atoms.empty()) {
        return {0.0, 0.0, 0.0};
    }
    std::array<double, 3> minv{std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::infinity()};
    for (const auto& atom : molecule.atoms) {
        minv[0] = std::min(minv[0], atom.x);
        minv[1] = std::min(minv[1], atom.y);
        minv[2] = std::min(minv[2], atom.z);
    }
    return {3.0 - minv[0], 3.0 - minv[1], 3.0 - minv[2]};
}

}  // namespace xpongecpp
