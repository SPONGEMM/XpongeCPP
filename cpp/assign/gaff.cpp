#include "core.hpp"
#include "amber_alternating.hpp"

#include <array>
#include <functional>
#include <stdexcept>
#include <unordered_set>

namespace xpongecpp {
namespace {

using RuleFunction = std::function<bool(std::uint32_t, const Assign&)>;

struct Rule {
    std::string name;
    RuleFunction match;
};

const std::vector<std::string>& gaff_rule_names() {
    static const std::vector<std::string> names{
        "cx", "cy", "c3", "c",  "cz", "cq", "cp", "ca", "cd", "cc", "cf", "ce", "cu", "cv", "c2",
        "cg", "c1", "hn", "ho", "hs", "hp", "hx", "hw", "h3", "h2", "h1", "hc", "h5", "h4", "ha",
        "f",  "cl", "br", "i",  "o",  "oh", "op", "oq", "os", "ni", "nj", "n",  "nk", "nl", "n4",
        "no", "na", "nm", "nn", "nh", "np", "nq", "n3", "nb", "nd", "nc", "nf", "ne", "n1", "n2",
        "s",  "s2", "sh", "sp", "sq", "ss", "sx", "s4", "sy", "s6", "pd", "pb", "pc", "pf", "pe",
        "p2", "px", "p4", "p3", "py", "p5",
    };
    return names;
}

const std::vector<std::string>& gaff2_rule_names() {
    static const std::vector<std::string> names{
        "cx", "cy", "c5", "c6", "c3", "cs", "c",  "cz", "cq", "cp", "ca", "cd", "cc", "cf", "ce", "cu", "cv",
        "c2", "cg", "c1", "hn", "ho", "hs", "hp", "hx", "hw", "h3", "h2", "h1", "hc", "h5", "h4", "ha", "f",
        "cl", "br", "i",  "o",  "oh", "op", "oq", "os", "ni", "nj", "ns", "nt", "n",  "nk", "nl", "nx", "ny",
        "nz", "n+", "n4", "no", "na", "nm", "nn", "nu", "nv", "nh", "np", "nq", "n5", "n6", "n7", "n8", "n9",
        "n3", "nb", "nd", "nc", "nf", "ne", "n1", "n2", "s",  "s2", "sh", "sp", "sq", "ss", "sx", "s4", "sy",
        "s6", "pd", "pb", "pc", "pf", "pe", "p2", "px", "p4", "p3", "py", "p5",
    };
    return names;
}

bool in_set(const std::string& value, const std::unordered_set<std::string>& values) {
    return values.count(value) != 0;
}

bool is_xx(const std::string& element) {
    static const std::unordered_set<std::string> values{"C", "N", "O", "P", "S"};
    return in_set(element, values);
}

bool is_xe(const std::string& element) {
    static const std::unordered_set<std::string> values{"N", "O", "F", "Cl", "Br", "S", "I"};
    return in_set(element, values);
}

bool has_any_marker(const Assign& assign, std::uint32_t atom, const std::vector<std::string>& markers) {
    for (const auto& marker : markers) {
        if (assign.has_atom_marker(atom, marker)) {
            return true;
        }
    }
    return false;
}

bool has_double_like_marker(const Assign& assign, std::uint32_t atom) {
    return assign.has_atom_marker(atom, "db") ||
           assign.has_atom_marker(atom, "DL") ||
           assign.has_atom_marker(atom, "DLB");
}

std::uint32_t first_neighbor(const Assign& assign, std::uint32_t atom) {
    if (assign.bonds[atom].empty()) {
        throw std::runtime_error("GAFF rule expected a bonded atom");
    }
    return assign.bonds[atom].begin()->first;
}

bool single_double_name(std::uint32_t atom, const Assign& assign, const std::string& name1, const std::string& name2) {
    for (const auto& [bonded_atom, bond_order] : assign.bonds[atom]) {
        if ((assign.atom_types[bonded_atom] == name1 && bond_order == 2) ||
            (assign.atom_types[bonded_atom] == name2 && bond_order == 1)) {
            return true;
        }
    }
    return false;
}

bool has_single_bond_to_conjugated_partner(std::uint32_t atom, const Assign& assign,
                                           const std::vector<std::string>& first_masks) {
    for (const auto& [bonded_atom, bond_order] : assign.bonds[atom]) {
        if (bond_order != 1) {
            continue;
        }
        if (assign.atom_judge(bonded_atom, first_masks)) {
            return true;
        }
        if (assign.atom_judge(bonded_atom, std::vector<std::string>{"S3", "S4", "P3", "P4"}) &&
            assign.has_atom_marker(bonded_atom, "db")) {
            return true;
        }
    }
    return false;
}

bool n3_amide_like(std::uint32_t atom, const Assign& assign) {
    if (!assign.atom_judge(atom, "N3")) {
        return false;
    }
    for (const auto& [bonded_atom, order] : assign.bonds[atom]) {
        (void)order;
        if (assign.atom_judge(bonded_atom, "C3")) {
            for (const auto& [bonded_atom_bonded, order2] : assign.bonds[bonded_atom]) {
                (void)order2;
                if (assign.atom_judge(bonded_atom_bonded, "O1") || assign.atom_judge(bonded_atom_bonded, "S1")) {
                    return true;
                }
            }
        }
    }
    return false;
}

bool n3_conjugated_ring_neighbor(std::uint32_t atom, const Assign& assign) {
    for (const auto& [bonded_atom, order] : assign.bonds[atom]) {
        (void)order;
        if (assign.has_atom_marker(bonded_atom, "DB") &&
            assign.atom_judge(bonded_atom, std::vector<std::string>{"C3", "N2", "P2"})) {
            return true;
        }
        if (has_any_marker(assign, bonded_atom, {"AR1", "AR2", "AR3"}) && is_xx(assign.elements[bonded_atom])) {
            return true;
        }
    }
    return false;
}

bool n2_conjugated(std::uint32_t atom, const Assign& assign) {
    if (!(assign.atom_judge(atom, "N2") && assign.has_atom_marker(atom, "sb") &&
          assign.has_atom_marker(atom, "db"))) {
        return false;
    }
    for (const auto& [bonded_atom, bond_order] : assign.bonds[atom]) {
        if (bond_order == 1 &&
            (assign.atom_judge(bonded_atom, std::vector<std::string>{"C3", "C2", "N2", "P2"}) ||
             (assign.atom_judge(bonded_atom, std::vector<std::string>{"S3", "S4", "P3", "P4"}) &&
              assign.has_atom_marker(bonded_atom, "db")))) {
            return true;
        }
        if (assign.atom_judge(bonded_atom, std::vector<std::string>{"C3", "C2", "N2", "P2"})) {
            for (const auto& [bonded_atom_bonded, order2] : assign.bonds[bonded_atom]) {
                (void)order2;
                if (assign.atom_judge(bonded_atom_bonded, std::vector<std::string>{"C3", "C2", "N2", "P2"})) {
                    return true;
                }
            }
        }
    }
    return false;
}

bool p2_sp2_conjugated(std::uint32_t atom, const Assign& assign, bool require_aromatic) {
    if (!(assign.atom_judge(atom, "P2") && assign.has_atom_marker(atom, "sb") &&
          assign.has_atom_marker(atom, "db"))) {
        return false;
    }
    if (require_aromatic && !has_any_marker(assign, atom, {"AR2", "AR3"})) {
        return false;
    }
    return has_single_bond_to_conjugated_partner(atom, assign, {"C3", "N2", "P2"});
}

bool p_or_s_db_conjugated(std::uint32_t atom, const Assign& assign, const std::string& mask) {
    if (!(assign.atom_judge(atom, mask) && assign.has_atom_marker(atom, "db"))) {
        return false;
    }
    return has_single_bond_to_conjugated_partner(atom, assign, {"C3", "N2", "P2"});
}

int attached_hydrogens(std::uint32_t atom, const Assign& assign) {
    int count = 0;
    for (const auto& [bonded_atom, bond_order] : assign.bonds[atom]) {
        (void)bond_order;
        if (assign.elements[bonded_atom] == "H") {
            ++count;
        }
    }
    return count;
}

bool has_amide_like_carbon_neighbor(std::uint32_t atom, const Assign& assign) {
    for (const auto& [bonded_atom, order] : assign.bonds[atom]) {
        (void)order;
        if (!assign.atom_judge(bonded_atom, "C3")) {
            continue;
        }
        for (const auto& [bonded_atom_bonded, order2] : assign.bonds[bonded_atom]) {
            (void)order2;
            if (assign.atom_judge(bonded_atom_bonded, "O1") || assign.atom_judge(bonded_atom_bonded, "S1")) {
                return true;
            }
        }
    }
    return false;
}

bool has_conjugated_amine_neighbor(std::uint32_t atom, const Assign& assign) {
    for (const auto& [bonded_atom, order] : assign.bonds[atom]) {
        (void)order;
        if (assign.has_atom_marker(bonded_atom, "DB") &&
            assign.atom_judge(bonded_atom, std::vector<std::string>{"C3", "N2", "P2"})) {
            return true;
        }
        if (has_any_marker(assign, bonded_atom, {"AR1", "AR2", "AR3"}) && is_xx(assign.elements[bonded_atom])) {
            return true;
        }
    }
    return false;
}

void assign_gaff_base_types(Assign& assign, const std::vector<Rule>& rules) {
    for (std::uint32_t atom = 0; atom < assign.elements.size(); ++atom) {
        bool found = false;
        for (const auto& current_rule : rules) {
            if (current_rule.match(atom, assign)) {
                assign.atom_types[atom] = current_rule.name;
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::runtime_error("No atom type found for assignment " + assign.name + " of atom #" + std::to_string(atom));
        }
    }
    apply_amber_alternating_type_adjustment(assign);
}

void apply_gaff2_refinements(Assign& assign) {
    for (std::uint32_t atom = 0; atom < assign.atom_count(); ++atom) {
        auto& atom_type = assign.atom_types[atom];
        if (atom_type == "c3") {
            if (assign.has_atom_marker(atom, "RG5")) {
                atom_type = "c5";
            } else if (assign.has_atom_marker(atom, "RG6")) {
                atom_type = "c6";
            }
            continue;
        }
        if (atom_type == "c") {
            for (const auto& [bonded_atom, order] : assign.bonds[atom]) {
                (void)order;
                if (assign.atom_judge(bonded_atom, "S1")) {
                    atom_type = "cs";
                    break;
                }
            }
            continue;
        }
        if (atom_type == "n4") {
            const auto hydrogens = attached_hydrogens(atom, assign);
            if (hydrogens == 1) {
                atom_type = "nx";
            } else if (hydrogens == 2) {
                atom_type = "ny";
            } else if (hydrogens == 3) {
                atom_type = "nz";
            } else if (hydrogens == 4) {
                atom_type = "n+";
            }
            continue;
        }
        if (atom_type == "n") {
            const auto hydrogens = attached_hydrogens(atom, assign);
            if (hydrogens == 1 && has_amide_like_carbon_neighbor(atom, assign)) {
                atom_type = "ns";
            } else if (hydrogens == 2 && has_amide_like_carbon_neighbor(atom, assign)) {
                atom_type = "nt";
            }
            continue;
        }
        if (atom_type == "nh") {
            const auto hydrogens = attached_hydrogens(atom, assign);
            if (hydrogens == 1 && has_conjugated_amine_neighbor(atom, assign)) {
                atom_type = "nu";
            } else if (hydrogens == 2 && has_conjugated_amine_neighbor(atom, assign)) {
                atom_type = "nv";
            }
            continue;
        }
        if (atom_type == "np" && attached_hydrogens(atom, assign) == 1) {
            atom_type = "n5";
            continue;
        }
        if (atom_type == "nq" && attached_hydrogens(atom, assign) == 1) {
            atom_type = "n6";
            continue;
        }
        if (atom_type == "n3") {
            const auto hydrogens = attached_hydrogens(atom, assign);
            if (hydrogens == 1) {
                atom_type = "n7";
            } else if (hydrogens == 2) {
                atom_type = "n8";
            } else if (hydrogens == 3) {
                atom_type = "n9";
            }
        }
    }
}

const std::vector<Rule>& gaff_rules() {
    static const std::vector<Rule> rules{
        {"cx", [](auto i, const auto& a) { return a.atom_judge(i, "C4") && a.has_atom_marker(i, "RG3"); }},
        {"cy", [](auto i, const auto& a) { return a.atom_judge(i, "C4") && a.has_atom_marker(i, "RG4"); }},
        {"c3", [](auto i, const auto& a) { return a.atom_judge(i, "C4"); }},
        {"c", [](auto i, const auto& a) {
             if (!a.atom_judge(i, "C3") || a.has_atom_marker(i, "AR1")) return false;
             for (const auto& [j, order] : a.bonds[i]) {
                 (void)order;
                 if (a.atom_judge(j, std::vector<std::string>{"O1", "S1"})) return true;
             }
             return false;
         }},
        {"cz", [](auto i, const auto& a) {
             if (!a.atom_judge(i, "C3")) return false;
             for (const auto& [j, order] : a.bonds[i]) {
                 (void)order;
                 if (!a.atom_judge(j, "N3")) return false;
             }
             return true;
         }},
        {"cq", [](auto i, const auto& a) { (void)i; (void)a; return false; }},
        {"cp", [](auto i, const auto& a) {
             bool tofind = a.atom_judge(i, "C3") && a.has_atom_marker(i, "AR1") && a.atom_marker_count(i, "RG6") == 1;
             for (const auto& [j, order] : a.bonds[i]) {
                 (void)order;
                 if (!tofind) break;
                 if (!is_xx(a.elements[j]) || !a.has_atom_marker(j, "AR1")) tofind = false;
             }
             return tofind;
         }},
        {"ca", [](auto i, const auto& a) { return a.atom_judge(i, "C3") && a.has_atom_marker(i, "AR1"); }},
        {"cd", [](auto i, const auto& a) { (void)i; (void)a; return false; }},
        {"cc", [](auto i, const auto& a) {
             return a.atom_judge(i, "C3") && a.has_atom_marker(i, "sb") && a.has_atom_marker(i, "db") &&
                    has_any_marker(a, i, {"AR2", "AR3"});
         }},
        {"cf", [](auto i, const auto& a) { (void)i; (void)a; return false; }},
        {"ce", [](auto i, const auto& a) {
             return a.atom_judge(i, "C3") && a.has_atom_marker(i, "sb") && a.has_atom_marker(i, "db") &&
                    has_single_bond_to_conjugated_partner(i, a, {"C3", "C2", "N2", "P2"});
         }},
        {"cu", [](auto i, const auto& a) { return a.atom_judge(i, "C3") && a.has_atom_marker(i, "RG3"); }},
        {"cv", [](auto i, const auto& a) { return a.atom_judge(i, "C3") && a.has_atom_marker(i, "RG4"); }},
        {"c2", [](auto i, const auto& a) { return a.atom_judge(i, "C3"); }},
        {"cg", [](auto i, const auto& a) {
             return a.atom_judge(i, "C2") && a.has_atom_marker(i, "sb") && a.has_atom_marker(i, "tb") &&
                    has_single_bond_to_conjugated_partner(i, a, {"C3", "C2", "N2", "P2", "N1"});
         }},
        {"c1", [](auto i, const auto& a) { return a.atom_judge(i, "C2") || a.atom_judge(i, "C1"); }},
        {"hn", [](auto i, const auto& a) { return a.atom_judge(i, "H1") && a.elements[first_neighbor(a, i)] == "N"; }},
        {"ho", [](auto i, const auto& a) { return a.atom_judge(i, "H1") && a.elements[first_neighbor(a, i)] == "O"; }},
        {"hs", [](auto i, const auto& a) { return a.atom_judge(i, "H1") && a.elements[first_neighbor(a, i)] == "S"; }},
        {"hp", [](auto i, const auto& a) { return a.atom_judge(i, "H1") && a.elements[first_neighbor(a, i)] == "P"; }},
        {"hx", [](auto i, const auto& a) {
             if (!a.atom_judge(i, "H1")) return false;
             for (const auto& [j, order] : a.bonds[i]) {
                 (void)order;
                 if (a.elements[j] == "C") {
                     for (const auto& [k, order2] : a.bonds[j]) {
                         (void)order2;
                         if (a.atom_judge(k, "N4")) return true;
                     }
                 }
             }
             return false;
         }},
        {"hw", [](auto i, const auto& a) {
             if (!a.atom_judge(i, "H1")) return false;
             for (const auto& [j, order] : a.bonds[i]) {
                 (void)order;
                 if (a.elements[j] == "O") {
                     for (const auto& [k, order2] : a.bonds[j]) {
                         (void)order2;
                         if (a.atom_judge(k, "H1")) return true;
                     }
                 }
             }
             return false;
         }},
        {"h3", [](auto i, const auto& a) {
             int count = 0;
             for (const auto& [j, order] : a.bonds[i]) {
                 (void)order;
                 if (a.atom_judge(j, "C4")) {
                     for (const auto& [k, order2] : a.bonds[j]) {
                         (void)order2;
                         if (is_xe(a.elements[k])) ++count;
                     }
                 }
             }
             return a.atom_judge(i, "H1") && count == 3;
         }},
        {"h2", [](auto i, const auto& a) {
             int count = 0;
             for (const auto& [j, order] : a.bonds[i]) {
                 (void)order;
                 if (a.atom_judge(j, "C4")) {
                     for (const auto& [k, order2] : a.bonds[j]) {
                         (void)order2;
                         if (is_xe(a.elements[k])) ++count;
                     }
                 }
             }
             return a.atom_judge(i, "H1") && count == 2;
         }},
        {"h1", [](auto i, const auto& a) {
             int count = 0;
             for (const auto& [j, order] : a.bonds[i]) {
                 (void)order;
                 if (a.atom_judge(j, "C4")) {
                     for (const auto& [k, order2] : a.bonds[j]) {
                         (void)order2;
                         if (is_xe(a.elements[k])) ++count;
                     }
                 }
             }
             return a.atom_judge(i, "H1") && count == 1;
         }},
        {"hc", [](auto i, const auto& a) { return a.atom_judge(i, "H1") && a.atom_judge(first_neighbor(a, i), "C4"); }},
        {"h5", [](auto i, const auto& a) {
             int count = 0;
             for (const auto& [j, order] : a.bonds[i]) {
                 (void)order;
                 if (a.atom_judge(j, "C3")) {
                     for (const auto& [k, order2] : a.bonds[j]) {
                         (void)order2;
                         if (is_xe(a.elements[k])) ++count;
                     }
                 }
             }
             return a.atom_judge(i, "H1") && count == 2;
         }},
        {"h4", [](auto i, const auto& a) {
             int count = 0;
             for (const auto& [j, order] : a.bonds[i]) {
                 (void)order;
                 if (a.atom_judge(j, "C3")) {
                     for (const auto& [k, order2] : a.bonds[j]) {
                         (void)order2;
                         if (is_xe(a.elements[k])) ++count;
                     }
                 }
             }
             return a.atom_judge(i, "H1") && count == 1;
         }},
        {"ha", [](auto i, const auto& a) { return a.atom_judge(i, "H1"); }},
        {"f", [](auto i, const auto& a) { return a.elements[i] == "F"; }},
        {"cl", [](auto i, const auto& a) { return a.elements[i] == "Cl"; }},
        {"br", [](auto i, const auto& a) { return a.elements[i] == "Br"; }},
        {"i", [](auto i, const auto& a) { return a.elements[i] == "I"; }},
        {"o", [](auto i, const auto& a) { return a.atom_judge(i, "O1"); }},
        {"oh", [](auto i, const auto& a) {
             if (!(a.atom_judge(i, "O2") || a.atom_judge(i, "O3"))) return false;
             for (const auto& [j, order] : a.bonds[i]) {
                 (void)order;
                 if (a.atom_judge(j, "H1")) return true;
             }
             return false;
         }},
        {"op", [](auto i, const auto& a) { return a.atom_judge(i, "O2") && a.has_atom_marker(i, "RG3"); }},
        {"oq", [](auto i, const auto& a) { return a.atom_judge(i, "O2") && a.has_atom_marker(i, "RG4"); }},
        {"os", [](auto i, const auto& a) { return a.elements[i] == "O"; }},
        {"ni", [](auto i, const auto& a) { return a.atom_judge(i, "N3") && a.has_atom_marker(i, "RG3") && n3_amide_like(i, a); }},
        {"nj", [](auto i, const auto& a) { return a.atom_judge(i, "N3") && a.has_atom_marker(i, "RG4") && n3_amide_like(i, a); }},
        {"n", [](auto i, const auto& a) { return n3_amide_like(i, a); }},
        {"nk", [](auto i, const auto& a) { return a.atom_judge(i, "N4") && a.has_atom_marker(i, "RG3"); }},
        {"nl", [](auto i, const auto& a) { return a.atom_judge(i, "N4") && a.has_atom_marker(i, "RG4"); }},
        {"n4", [](auto i, const auto& a) { return a.atom_judge(i, "N4"); }},
        {"no", [](auto i, const auto& a) {
             int count = 0;
             if (a.atom_judge(i, "N3")) {
                 for (const auto& [j, order] : a.bonds[i]) {
                     (void)order;
                     if (a.atom_judge(j, "O1")) ++count;
                 }
             }
             return count == 2;
         }},
        {"na", [](auto i, const auto& a) {
             return a.atom_judge(i, "N3") && has_any_marker(a, i, {"AR1", "AR2", "AR3"});
         }},
        {"nm", [](auto i, const auto& a) {
             return a.atom_judge(i, "N3") && a.has_atom_marker(i, "RG3") && n3_conjugated_ring_neighbor(i, a);
         }},
        {"nn", [](auto i, const auto& a) {
             return a.atom_judge(i, "N3") && a.has_atom_marker(i, "RG4") && n3_conjugated_ring_neighbor(i, a);
         }},
        {"nh", [](auto i, const auto& a) { return a.atom_judge(i, "N3") && n3_conjugated_ring_neighbor(i, a); }},
        {"np", [](auto i, const auto& a) { return a.atom_judge(i, "N3") && a.has_atom_marker(i, "RG3"); }},
        {"nq", [](auto i, const auto& a) { return a.atom_judge(i, "N3") && a.has_atom_marker(i, "RG4"); }},
        {"n3", [](auto i, const auto& a) { return a.atom_judge(i, "N3"); }},
        {"nb", [](auto i, const auto& a) { return a.atom_judge(i, "N2") && a.has_atom_marker(i, "AR1"); }},
        {"nd", [](auto i, const auto& a) { (void)i; (void)a; return false; }},
        {"nc", [](auto i, const auto& a) { return n2_conjugated(i, a) && has_any_marker(a, i, {"AR2", "AR3"}); }},
        {"nf", [](auto i, const auto& a) { (void)i; (void)a; return false; }},
        {"ne", [](auto i, const auto& a) {
             return a.atom_judge(i, "N2") && a.has_atom_marker(i, "sb") && has_double_like_marker(a, i) &&
                    has_single_bond_to_conjugated_partner(i, a, {"C3", "C2", "N2", "P2"});
         }},
        {"n1", [](auto i, const auto& a) {
             return a.atom_judge(i, "N1") ||
                    (a.atom_judge(i, "N2") &&
                     ((a.has_atom_marker(i, "sb") && a.has_atom_marker(i, "tb")) || a.atom_marker_count(i, "db") == 2));
         }},
        {"n2", [](auto i, const auto& a) { return a.atom_judge(i, "N2"); }},
        {"s", [](auto i, const auto& a) { return a.atom_judge(i, "S1"); }},
        {"s2", [](auto i, const auto& a) {
             return a.atom_judge(i, "S2") && (a.has_atom_marker(i, "DB") || a.has_atom_marker(i, "TB"));
         }},
        {"sh", [](auto i, const auto& a) {
             if (!a.atom_judge(i, "S2")) return false;
             for (const auto& [j, order] : a.bonds[i]) {
                 (void)order;
                 if (a.atom_judge(j, "H1")) return true;
             }
             return false;
         }},
        {"sp", [](auto i, const auto& a) { return a.atom_judge(i, "S2") && a.has_atom_marker(i, "RG3"); }},
        {"sq", [](auto i, const auto& a) { return a.atom_judge(i, "S2") && a.has_atom_marker(i, "RG4"); }},
        {"ss", [](auto i, const auto& a) { return a.atom_judge(i, "S2"); }},
        {"sx", [](auto i, const auto& a) { return p_or_s_db_conjugated(i, a, "S3"); }},
        {"s4", [](auto i, const auto& a) { return a.atom_judge(i, "S3"); }},
        {"sy", [](auto i, const auto& a) { return p_or_s_db_conjugated(i, a, "S4"); }},
        {"s6", [](auto i, const auto& a) { return a.atom_judge(i, "S4") || a.atom_judge(i, "S5") || a.atom_judge(i, "S6"); }},
        {"pd", [](auto i, const auto& a) { (void)i; (void)a; return false; }},
        {"pb", [](auto i, const auto& a) { return a.atom_judge(i, "P2") && a.has_atom_marker(i, "AR1"); }},
        {"pc", [](auto i, const auto& a) { return p2_sp2_conjugated(i, a, true); }},
        {"pf", [](auto i, const auto& a) { (void)i; (void)a; return false; }},
        {"pe", [](auto i, const auto& a) { return p2_sp2_conjugated(i, a, false); }},
        {"p2", [](auto i, const auto& a) { return a.atom_judge(i, "P1") || a.atom_judge(i, "P2"); }},
        {"px", [](auto i, const auto& a) { return p_or_s_db_conjugated(i, a, "P3"); }},
        {"p4", [](auto i, const auto& a) {
             if (!(a.atom_judge(i, "P3") && a.has_atom_marker(i, "db"))) return false;
             for (const auto& [j, order] : a.bonds[i]) {
                 (void)order;
                 if (a.atom_judge(j, "O1") || a.atom_judge(j, "S1")) return true;
             }
             return false;
         }},
        {"p3", [](auto i, const auto& a) { return a.atom_judge(i, "P3"); }},
        {"py", [](auto i, const auto& a) { return p_or_s_db_conjugated(i, a, "P4"); }},
        {"p5", [](auto i, const auto& a) { return a.atom_judge(i, "P4") || a.atom_judge(i, "P5") || a.atom_judge(i, "P6"); }},
    };
    return rules;
}

}  // namespace

std::vector<std::string> implemented_gaff_assign_types() {
    return gaff_rule_names();
}

std::vector<std::string> implemented_gaff2_assign_types() {
    return gaff2_rule_names();
}

void Assign::determine_atom_type(const std::string& rule) {
    if (rule == "sybyl" || rule == "SYBYL") {
        for (std::size_t i = 0; i < elements.size(); ++i) {
            atom_types[i] = elements[i] + element_details[i];
        }
        return;
    }
    const bool is_gaff = rule == "gaff" || rule == "GAFF";
    const bool is_gaff2 = rule == "gaff2" || rule == "GAFF2";
    if (!is_gaff && !is_gaff2) {
        for (std::size_t i = 0; i < elements.size(); ++i) {
            atom_types[i] = elements[i];
        }
        return;
    }
    if (!built) {
        determine_ring_and_bond_type();
    }
    const auto& rules = gaff_rules();
    assign_gaff_base_types(*this, rules);
    if (is_gaff2) {
        apply_gaff2_refinements(*this);
    }
}

}  // namespace xpongecpp
