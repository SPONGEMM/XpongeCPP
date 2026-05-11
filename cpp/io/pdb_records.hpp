#pragma once

#include <string>

namespace xpongecpp {

std::string pdb_trim_copy(const std::string& input);
int pdb_hy36_decode(int width, const std::string& field);
std::string pdb_hy36_field(int width, int value);

}  // namespace xpongecpp
