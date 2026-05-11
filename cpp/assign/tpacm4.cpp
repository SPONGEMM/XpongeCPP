#include "core.hpp"

#include <algorithm>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace xpongecpp {
namespace {

struct Tpacm4Table {
    std::unordered_map<std::string, double> charge_by_type;
    std::vector<std::pair<std::string, double>> ordered;
};

std::vector<std::string> split_ws(const std::string& text) {
    std::istringstream input(text);
    std::vector<std::string> words;
    std::string word;
    while (input >> word) {
        words.push_back(word);
    }
    return words;
}

Tpacm4Table parse_tpacm4_table(const std::string& atom_type_table, const std::string& charge_table) {
    const auto atom_types = split_ws(atom_type_table);
    const auto charges = split_ws(charge_table);
    if (atom_types.size() != charges.size()) {
        throw std::invalid_argument("TPACM4 atom type and charge table lengths differ");
    }
    Tpacm4Table table;
    table.ordered.reserve(atom_types.size());
    for (std::size_t i = 0; i < atom_types.size(); ++i) {
        const double charge = std::stod(charges[i]);
        table.charge_by_type[atom_types[i]] = charge;
        table.ordered.emplace_back(atom_types[i], charge);
    }
    return table;
}

std::string tpacm4_element_code(const std::string& element) {
    if (element == "Cl" || element == "Br") {
        return "B";
    }
    if (element == "C" || element == "N" || element == "O" || element == "F" || element == "P" ||
        element == "S" || element == "H") {
        return element;
    }
    throw std::invalid_argument("TPACM4 does not support element " + element);
}

std::vector<std::vector<std::string>> unique_permutations(std::vector<std::string> values) {
    std::vector<std::vector<std::string>> out;
    std::sort(values.begin(), values.end());
    do {
        out.push_back(values);
    } while (std::next_permutation(values.begin(), values.end()));
    return out;
}

std::string join_strings(const std::vector<std::string>& values) {
    std::string out;
    for (const auto& value : values) {
        out += value;
    }
    return out;
}

std::string tpacm4_suffix(const Assign& assign, std::uint32_t atom) {
    if (assign.elements[atom] == "N") {
        int order_sum = 0;
        for (const auto& [neighbor, order] : assign.bonds[atom]) {
            (void)neighbor;
            order_sum += order;
        }
        if (order_sum == 4) {
            return "NH4";
        }
    }
    if (assign.elements[atom] == "O") {
        for (const auto& [neighbor, order] : assign.bonds[atom]) {
            if (order != 1 || !assign.atom_judge(neighbor, "C3")) {
                continue;
            }
            for (const auto& [second, second_order] : assign.bonds[neighbor]) {
                if (assign.elements[second] == "O" && second_order == 2) {
                    return "EST";
                }
            }
        }
    }
    if (assign.elements[atom] == "N") {
        for (const auto& [neighbor, order] : assign.bonds[atom]) {
            if (order != 1 || !assign.atom_judge(neighbor, "C3")) {
                continue;
            }
            for (const auto& [second, second_order] : assign.bonds[neighbor]) {
                if (assign.elements[second] == "O" && second_order == 2) {
                    return "ACM";
                }
            }
        }
    }
    return "";
}

int tpacm4_type_distance(const std::string& s1, const std::string& s2) {
    if (s1.size() != s2.size()) {
        return 0;
    }
    int out = 0;
    for (std::size_t i = 0; i < s1.size(); ++i) {
        if (i % 3 != 2 && s1[i] != s2[i]) {
            return 0;
        }
        if (i % 3 == 2 && s1[i] != s2[i]) {
            ++out;
        }
    }
    return out;
}

std::unordered_map<std::string, int> extra_string_order() {
    const std::vector<std::string> names{"NH4", "3MR", "4MR", "5MR", "OEW", "CC4", "OED",
                                         "CO2", "CN2", "CN3", "OCO", "oPY", "mPY", "PPY",
                                         "oFU", "mFU", "oTF", "mTF", "OXX"};
    std::unordered_map<std::string, int> order;
    for (std::size_t i = 0; i < names.size(); ++i) {
        order[names[i]] = static_cast<int>(i);
    }
    return order;
}

void ring_process(const Assign& assign,
                  const std::vector<std::uint32_t>& ring,
                  std::vector<std::vector<std::string>>& extra_strings) {
    if (ring.size() < 6) {
        for (const auto atom : ring) {
            extra_strings[atom].push_back(std::to_string(ring.size()) + "MR");
        }
    }
    const bool aromatic = ring.size() <= 6 &&
                          std::all_of(ring.begin(), ring.end(), [&](std::uint32_t atom) {
                              return assign.has_atom_marker(atom, "AR0");
                          });
    if (!aromatic) {
        return;
    }
    for (std::size_t i = 0; i < ring.size(); ++i) {
        const auto atom = ring[i];
        std::string ring_name;
        if (assign.elements[atom] == "N") {
            ring_name = "PY";
        } else if (assign.elements[atom] == "O") {
            ring_name = "FU";
        } else if (assign.elements[atom] == "S") {
            ring_name = "TF";
        }
        if (ring_name.empty()) {
            continue;
        }
        for (const auto offset : {ring.size() - 1, std::size_t{1}}) {
            const auto btom = ring[(i + offset) % ring.size()];
            for (const auto& [ctom, order] : assign.bonds[btom]) {
                (void)order;
                if (!assign.has_bond_marker(btom, ctom, "ar")) {
                    extra_strings[ctom].push_back("o" + ring_name);
                }
            }
        }
        for (const auto offset : {ring.size() - 2, std::size_t{2}}) {
            const auto btom = ring[(i + offset) % ring.size()];
            for (const auto& [ctom, order] : assign.bonds[btom]) {
                (void)order;
                if (!assign.has_bond_marker(btom, ctom, "ar")) {
                    extra_strings[ctom].push_back("m" + ring_name);
                }
            }
        }
        if (ring.size() == 6) {
            extra_strings[ring[(i + ring.size() - 3) % ring.size()]].push_back("P" + ring_name);
        }
    }
}

std::vector<std::string> find_extra_strings(const std::vector<std::string>& atom_type_alls, const Assign& assign) {
    std::vector<std::vector<std::string>> extra_strings(assign.atom_count());
    for (const auto& ring : assign.rings) {
        ring_process(assign, ring, extra_strings);
    }
    for (std::uint32_t i = 0; i < assign.atom_count(); ++i) {
        if (std::find(extra_strings[i].begin(), extra_strings[i].end(), "3MR") != extra_strings[i].end() ||
            std::find(extra_strings[i].begin(), extra_strings[i].end(), "4MR") != extra_strings[i].end() ||
            std::find(extra_strings[i].begin(), extra_strings[i].end(), "5MR") != extra_strings[i].end()) {
            continue;
        }
        const auto& type = atom_type_alls[i];
        if (type.find("NH4") != std::string::npos || type.find("NO2") != std::string::npos ||
            type.find("CN3") != std::string::npos || type.find("SO2") != std::string::npos) {
            for (const auto& [j, order] : assign.bonds[i]) {
                (void)order;
                if (assign.elements[j] != "C") {
                    continue;
                }
                for (const auto& [k, korder] : assign.bonds[j]) {
                    (void)korder;
                    if (k != i && assign.elements[k] == "C") {
                        extra_strings[k].push_back("OEW");
                    }
                }
            }
        } else if (((type.find("OH1") != std::string::npos || type.find("OC1") != std::string::npos) &&
                    assign.atom_judge(i, "O2")) ||
                   ((type.find("NH1") != std::string::npos || type.find("NC1") != std::string::npos) &&
                    assign.atom_judge(i, "N3"))) {
            for (const auto& [j, order] : assign.bonds[i]) {
                (void)order;
                if (assign.elements[j] != "C" || atom_type_alls[j].find("CO2") != std::string::npos) {
                    continue;
                }
                for (const auto& [k, korder] : assign.bonds[j]) {
                    (void)korder;
                    if (k != i && assign.elements[k] == "C") {
                        extra_strings[k].push_back("OED");
                    }
                }
            }
        } else if (assign.elements[i] == "F") {
            for (const auto& [j, order] : assign.bonds[i]) {
                (void)order;
                if (assign.elements[j] != "C") {
                    continue;
                }
                for (const auto& [k, korder] : assign.bonds[j]) {
                    (void)korder;
                    if (k != i && assign.elements[k] == "C") {
                        extra_strings[k].push_back("OXX");
                    }
                }
            }
        } else if (type.find("SH1") != std::string::npos || type.find("SC1") != std::string::npos) {
            for (const auto& [j, order] : assign.bonds[i]) {
                (void)order;
                if (assign.elements[j] != "C") {
                    continue;
                }
                for (const auto& [k, korder] : assign.bonds[j]) {
                    (void)korder;
                    if (k != i && assign.elements[k] == "C") {
                        extra_strings[k].push_back("OSH");
                    }
                }
            }
        }
    }

    const auto sort_order = extra_string_order();
    std::vector<std::string> out(assign.atom_count());
    for (std::uint32_t i = 0; i < assign.atom_count(); ++i) {
        if (assign.elements[i] == "H" && !assign.bonds[i].empty()) {
            const auto j = assign.bonds[i].begin()->first;
            if (std::find(extra_strings[j].begin(), extra_strings[j].end(), "OED") != extra_strings[j].end()) {
                extra_strings[i].push_back("OED");
            }
            if (std::find(extra_strings[j].begin(), extra_strings[j].end(), "OEW") != extra_strings[j].end()) {
                extra_strings[i].push_back("OEW");
            }
            for (const auto& [k, order] : assign.bonds[j]) {
                (void)order;
                if (atom_type_alls[k].find("CO2") != std::string::npos) {
                    extra_strings[i].push_back("CO2");
                } else if (atom_type_alls[k].find("CN2") != std::string::npos) {
                    extra_strings[i].push_back("CN2");
                } else if (atom_type_alls[k].find("NC3") != std::string::npos ||
                           atom_type_alls[k].find("CN3") != std::string::npos) {
                    extra_strings[i].push_back("CN3");
                } else if (atom_type_alls[k].find("NH4") != std::string::npos) {
                    extra_strings[i].push_back("NH4");
                }
            }
        }
        std::sort(extra_strings[i].begin(), extra_strings[i].end(), [&](const auto& left, const auto& right) {
            const int left_order = sort_order.count(left) == 0 ? 999 : sort_order.at(left);
            const int right_order = sort_order.count(right) == 0 ? 999 : sort_order.at(right);
            if (left_order != right_order) {
                return left_order < right_order;
            }
            return left < right;
        });
        out[i] = join_strings(extra_strings[i]);
    }
    return out;
}

}  // namespace

void Assign::calculate_tpacm4_charge(const std::string& atom_type_table,
                                     const std::string& charge_table,
                                     int total_charge) {
    const auto table = parse_tpacm4_table(atom_type_table, charge_table);
    kekulize();
    std::vector<std::string> suffixes(atom_count());
    for (std::uint32_t atom = 0; atom < atom_count(); ++atom) {
        suffixes[atom] = tpacm4_suffix(*this, atom);
    }

    std::vector<std::string> atom_type_alls;
    atom_type_alls.reserve(atom_count());
    std::vector<double> new_charges;
    new_charges.reserve(atom_count());
    for (std::uint32_t atom = 0; atom < atom_count(); ++atom) {
        std::vector<std::string> atom_type;
        for (const auto& [neighbor, order_in] : bonds[atom]) {
            int order = order_in;
            if (has_bond_marker(atom, neighbor, "ar")) {
                order = 4;
            }
            atom_type.push_back(tpacm4_element_code(elements[atom]) + tpacm4_element_code(elements[neighbor]) +
                                std::to_string(order));
        }

        bool found = false;
        for (const auto& one_pos : unique_permutations(atom_type)) {
            const auto temp_atom_type = join_strings(one_pos) + suffixes[atom];
            const auto type_it = table.charge_by_type.find(temp_atom_type);
            if (type_it != table.charge_by_type.end()) {
                atom_type_alls.push_back(temp_atom_type);
                new_charges.push_back(type_it->second);
                found = true;
                break;
            }
        }
        const auto total_length = atom_type.size();
        for (std::size_t dis = 0; dis < total_length && !found; ++dis) {
            for (const auto& one_pos : unique_permutations(atom_type)) {
                if (found) {
                    break;
                }
                const auto temp_atom_type = join_strings(one_pos);
                for (const auto& [type, charge] : table.ordered) {
                    if (tpacm4_type_distance(type, temp_atom_type) == static_cast<int>(dis)) {
                        atom_type_alls.push_back(type);
                        new_charges.push_back(charge);
                        found = true;
                        break;
                    }
                }
            }
        }
        if (!found) {
            atom_type_alls.push_back("XXX");
            new_charges.push_back(0.0);
        }
    }

    const auto extra_strings = find_extra_strings(atom_type_alls, *this);
    for (std::uint32_t atom = 0; atom < atom_count(); ++atom) {
        std::vector<std::string> neighbors;
        neighbors.reserve(bonds[atom].size());
        for (const auto& [neighbor, order] : bonds[atom]) {
            (void)order;
            neighbors.push_back(atom_type_alls[neighbor]);
        }
        std::sort(neighbors.begin(), neighbors.end(), [](std::string left, std::string right) {
            std::replace(left.begin(), left.end(), 'H', 'Z');
            std::replace(right.begin(), right.end(), 'H', 'Z');
            return left < right;
        });
        for (const auto& one_pos : unique_permutations(neighbors)) {
            const auto temp_atom_type = atom_type_alls[atom] + join_strings(one_pos) + extra_strings[atom];
            const auto type_it = table.charge_by_type.find(temp_atom_type);
            if (type_it != table.charge_by_type.end()) {
                new_charges[atom] = type_it->second;
                break;
            }
        }
    }

    const double sum = std::accumulate(new_charges.begin(), new_charges.end(), 0.0);
    const double delta = (sum - static_cast<double>(total_charge)) / static_cast<double>(atom_count());
    for (auto& charge : new_charges) {
        charge -= delta;
    }
    charges = std::move(new_charges);
}

}  // namespace xpongecpp
