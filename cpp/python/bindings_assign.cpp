#include "bindings_internal.hpp"

namespace xpongecpp {
namespace {

std::shared_ptr<Assign> get_assignment_from_mol2_object(py::object source, py::object total_charge) {
    std::optional<int> charge;
    bool charge_from_sum = false;
    if (!total_charge.is_none()) {
        if (py::isinstance<py::str>(total_charge)) {
            const auto value = py::cast<std::string>(total_charge);
            if (value == "sum") {
                charge_from_sum = true;
            } else {
                throw std::invalid_argument("total_charge string must be 'sum'");
            }
        } else {
            charge = py::cast<int>(total_charge);
        }
    }
    return std::make_shared<Assign>(get_assignment_from_mol2_text(read_python_input(source), charge, charge_from_sum));
}

std::shared_ptr<Assign> get_assignment_from_xyz_object(py::object source) {
    return std::make_shared<Assign>(get_assignment_from_xyz_text(read_python_input(source)));
}

std::shared_ptr<Assign> get_assignment_from_pdb_object(py::object source) {
    return std::make_shared<Assign>(get_assignment_from_pdb_text(read_python_input(source)));
}

void save_assignment_mol2_object(const Assign& assignment, const std::string& filename,
                                 const std::string& residue_name) {
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("failed to open assignment MOL2 output: " + filename);
    }
    out << assignment_to_mol2_text(assignment, residue_name);
}

void save_assignment_pdb_object(const Assign& assignment, const std::string& filename,
                                const std::string& residue_name) {
    std::ofstream out(filename);
    if (!out) {
        throw std::runtime_error("failed to open assignment PDB output: " + filename);
    }
    out << assignment_to_pdb_text(assignment, residue_name);
}

}  // namespace

void bind_assign_module(py::module_& m) {
    py::class_<Assign, std::shared_ptr<Assign>>(m, "Assign")
        .def(py::init<std::string>(), py::arg("name") = "ASN")
        .def_readwrite("name", &Assign::name)
        .def_property_readonly("atom_count", &Assign::atom_count)
        .def_property_readonly("bond_count", &Assign::bond_count)
        .def_property_readonly("atoms", [](const Assign& self) { return self.elements; })
        .def_readwrite("built", &Assign::built)
        .def_readonly("element_details", &Assign::element_details)
        .def_readonly("names", &Assign::names)
        .def_readonly("coordinates", &Assign::coordinates)
        .def_readonly("charges", &Assign::charges)
        .def_readonly("formal_charges", &Assign::formal_charges)
        .def_readonly("bonds", &Assign::bonds)
        .def_readonly("bond_sequence", &Assign::bond_sequence)
        .def_readonly("bond_markers", &Assign::bond_markers)
        .def_readonly("atom_markers", &Assign::atom_markers)
        .def_readonly("atom_types", &Assign::atom_types)
        .def("add_atom", &Assign::add_atom, py::arg("element"), py::arg("x"), py::arg("y"), py::arg("z"),
             py::arg("name") = "", py::arg("charge") = 0.0)
        .def("Add_Atom", &Assign::add_atom, py::arg("element"), py::arg("x"), py::arg("y"), py::arg("z"),
             py::arg("name") = "", py::arg("charge") = 0.0)
        .def("add_bond", &Assign::add_bond, py::arg("atom1"), py::arg("atom2"), py::arg("order") = 1)
        .def("delete_bond", &Assign::delete_bond, py::arg("atom1"), py::arg("atom2"))
        .def("delete_atom", &Assign::delete_atom, py::arg("atom"))
        .def("set_charge", &Assign::set_charge, py::arg("atom"), py::arg("charge"))
        .def("set_charges", &Assign::set_charges, py::arg("charges"))
        .def("set_formal_charge", &Assign::set_formal_charge, py::arg("atom"), py::arg("charge"))
        .def("set_coordinate", &Assign::set_coordinate, py::arg("atom"), py::arg("x"), py::arg("y"), py::arg("z"))
        .def("set_atom_type", &Assign::set_atom_type, py::arg("atom"), py::arg("atom_type"))
        .def("atom_judge", py::overload_cast<std::uint32_t, const std::string&>(&Assign::atom_judge, py::const_),
             py::arg("atom"), py::arg("mask"))
        .def("atom_judge",
             py::overload_cast<std::uint32_t, const std::vector<std::string>&>(&Assign::atom_judge, py::const_),
             py::arg("atom"), py::arg("masks"))
        .def("add_atom_marker", &Assign::add_atom_marker, py::arg("atom"), py::arg("marker"))
        .def("has_atom_marker", &Assign::has_atom_marker, py::arg("atom"), py::arg("marker"))
        .def("atom_marker_count", &Assign::atom_marker_count, py::arg("atom"), py::arg("marker"))
        .def("has_bond_marker", &Assign::has_bond_marker, py::arg("atom1"), py::arg("atom2"), py::arg("marker"))
        .def("add_bond_marker", &Assign::add_bond_marker, py::arg("atom1"), py::arg("atom2"),
             py::arg("marker"), py::arg("only1") = false)
        .def("determine_connectivity", &Assign::determine_connectivity, py::arg("simple_cutoff"))
        .def("determine_bond_order", &Assign::determine_bond_order, py::arg("check_formal_charge") = true,
             py::arg("total_charge") = py::none())
        .def("_determine_bond_order_custom",
             [](Assign& self,
                bool check_formal_charge,
                py::object total_charge_object,
                int max_step,
                int max_stat,
                const std::vector<std::vector<std::pair<int, int>>>& penalty_scores,
                py::object extra_criteria) {
                 std::optional<int> total_charge;
                 if (!total_charge_object.is_none()) {
                     total_charge = py::cast<int>(total_charge_object);
                 }
                 std::function<bool(const Assign&)> criteria;
                 if (!extra_criteria.is_none()) {
                     criteria = [&self, extra_criteria](const Assign&) {
                         return py::cast<bool>(extra_criteria(py::cast(&self, py::return_value_policy::reference)));
                     };
                 }
                 return self.determine_bond_order_custom(check_formal_charge, total_charge, max_step, max_stat,
                                                         penalty_scores, criteria);
             },
             py::arg("check_formal_charge"), py::arg("total_charge"), py::arg("max_step"), py::arg("max_stat"),
             py::arg("penalty_scores"), py::arg("extra_criteria") = py::none())
        .def("determine_ring_and_bond_type", &Assign::determine_ring_and_bond_type)
        .def("kekulize", &Assign::kekulize)
        .def("determine_atom_type", &Assign::determine_atom_type, py::arg("rule"))
        .def("_calculate_tpacm4", &Assign::calculate_tpacm4_charge, py::arg("atom_type_table"),
             py::arg("charge_table"), py::arg("total_charge"))
        .def("to_residuetype", &Assign::to_residuetype, py::arg("name"))
        .def("to_molecule", [](const std::shared_ptr<Assign>& self, const std::string& name) {
            return std::make_shared<Molecule>(self->to_molecule(name));
        }, py::arg("name"))
        .def("save_as_mol2", &save_assignment_mol2_object, py::arg("filename"), py::arg("residue_name") = "MOL")
        .def("Save_As_Mol2", &save_assignment_mol2_object, py::arg("filename"), py::arg("residue_name") = "MOL")
        .def("save_as_pdb", &save_assignment_pdb_object, py::arg("filename"), py::arg("residue_name") = "MOL")
        .def("Save_As_PDB", &save_assignment_pdb_object, py::arg("filename"), py::arg("residue_name") = "MOL");

    m.def("get_assignment_from_mol2", &get_assignment_from_mol2_object, py::arg("source"),
          py::arg("total_charge") = py::none());
    m.def("get_assignment_from_xyz", &get_assignment_from_xyz_object, py::arg("source"));
    m.def("get_assignment_from_pdb", &get_assignment_from_pdb_object, py::arg("source"));
    m.def("get_assignment_from_residuetype",
          [](const ResidueType& residue_type) {
              return std::make_shared<Assign>(get_assignment_from_residuetype(residue_type));
          }, py::arg("residue_type"));
    m.def("implemented_gaff_assign_types", &implemented_gaff_assign_types);
    m.def("implemented_gaff2_assign_types", &implemented_gaff2_assign_types);
    m.def("generate_resp_mk_grid", &generate_resp_mk_grid,
          py::arg("atoms"), py::arg("atom_coordinates_bohr"), py::arg("area_density") = 1.0,
          py::arg("layer") = 4, py::arg("radius") = std::unordered_map<std::string, double>{});
    m.def("fit_resp_from_esp_cpp", &fit_resp_from_esp_cpp,
          py::arg("assign"), py::arg("atom_coordinates_bohr"), py::arg("nuclear_charges"),
          py::arg("grid_points_bohr"), py::arg("esp_values_au"), py::arg("charge"),
          py::arg("extra_equivalence") = std::vector<std::vector<int>>{},
          py::arg("a1") = 0.0005, py::arg("a2") = 0.001,
          py::arg("two_stage") = true, py::arg("only_esp") = false);
    m.def("fit_resp_from_esp_cpp_debug",
          [](const Assign& assign,
             const std::vector<std::array<double, 3>>& atom_coordinates_bohr,
             const std::vector<double>& nuclear_charges,
             const std::vector<std::array<double, 3>>& grid_points_bohr,
             const std::vector<double>& esp_values_au,
             int charge,
             const std::vector<std::vector<int>>& extra_equivalence,
             double a1,
             double a2,
             bool two_stage,
             bool only_esp) {
              const auto result = fit_resp_from_esp_cpp_debug(
                  assign,
                  atom_coordinates_bohr,
                  nuclear_charges,
                  grid_points_bohr,
                  esp_values_au,
                  charge,
                  extra_equivalence,
                  a1,
                  a2,
                  two_stage,
                  only_esp
              );
              py::dict timings;
              for (const auto& item : result.timings) {
                  timings[py::str(item.first)] = item.second;
              }
              py::dict out;
              out["esp_charges"] = result.esp_charges;
              out["stage1_charges"] = result.stage1_charges;
              out["final_charges"] = result.final_charges;
              out["stage2_restrained_groups"] = result.stage2_restrained_groups;
              out["timings"] = timings;
              return out;
          },
          py::arg("assign"), py::arg("atom_coordinates_bohr"), py::arg("nuclear_charges"),
          py::arg("grid_points_bohr"), py::arg("esp_values_au"), py::arg("charge"),
          py::arg("extra_equivalence") = std::vector<std::vector<int>>{},
          py::arg("a1") = 0.0005, py::arg("a2") = 0.001,
          py::arg("two_stage") = true, py::arg("only_esp") = false);
}

}  // namespace xpongecpp
