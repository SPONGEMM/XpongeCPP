#include "pdb_records.hpp"

#include <cctype>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace xpongecpp {

std::string pdb_trim_copy(const std::string& input) {
    const auto first = input.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return "";
    }
    const auto last = input.find_last_not_of(" \t\r\n");
    return input.substr(first, last - first + 1);
}

int pdb_hy36_decode(int width, const std::string& field) {
    const auto text = pdb_trim_copy(field);
    if (text.empty()) {
        return 0;
    }
    if (text[0] == '-' || text[0] == '+' || std::isdigit(static_cast<unsigned char>(text[0]))) {
        return std::stoi(text);
    }
    const std::string digits_upper = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    const std::string digits_lower = "0123456789abcdefghijklmnopqrstuvwxyz";
    const bool lower = std::islower(static_cast<unsigned char>(field[0]));
    const auto& digits = lower ? digits_lower : digits_upper;
    int value = 0;
    for (const char c : field) {
        const auto pos = digits.find(c);
        if (pos == std::string::npos) {
            throw std::invalid_argument("invalid hybrid-36 number");
        }
        value = value * 36 + static_cast<int>(pos);
    }
    int power = 1;
    for (int i = 1; i < width; ++i) {
        power *= 36;
    }
    if (lower) {
        return value + 16 * power + static_cast<int>(std::pow(10, width));
    }
    return value - 10 * power + static_cast<int>(std::pow(10, width));
}

std::string pdb_hy36_field(int width, int value) {
    const int decimal_limit = static_cast<int>(std::pow(10, width));
    if (value > -decimal_limit / 10 && value < decimal_limit) {
        std::ostringstream out;
        out << std::setw(width) << value;
        return out.str();
    }
    const std::string digits = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    int offset = value + 10 * static_cast<int>(std::pow(36, width - 1)) - decimal_limit;
    std::string out(width, '0');
    for (int i = width - 1; i >= 0; --i) {
        out[static_cast<std::size_t>(i)] = digits[static_cast<std::size_t>(offset % 36)];
        offset /= 36;
    }
    return out;
}

}  // namespace xpongecpp
