#pragma once

#include "core.hpp"

#include <array>

namespace xpongecpp {

std::array<double, 6> sponge_coordinate_box_for_export(const Molecule& molecule);
std::array<double, 3> sponge_coordinate_shift_for_export(const Molecule& molecule);

}  // namespace xpongecpp
