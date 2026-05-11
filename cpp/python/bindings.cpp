#include "core.hpp"

#include <fstream>
#include <memory>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace xpongecpp;

namespace {

std::string read_python_input(py::object source) {
    if (py::hasattr(source, "read")) {
        return py::cast<std::string>(source.attr("read")());
    }
    const auto path = py::cast<std::string>(source);
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("failed to open file: " + path);
    }
    std::ostringstream buffer;
    buffer << input.rdbuf();
    return buffer.str();
}

struct AtomView {
    std::shared_ptr<Molecule> molecule;
    AtomId id{0};

    const Atom& get() const { return molecule->atom(id); }
};

struct ResidueView {
    std::shared_ptr<Molecule> molecule;
    ResidueId id{0};

    const Residue& get() const { return molecule->residue(id); }

    AtomView atom(std::uint32_t local_index) const {
        const auto& residue = get();
        if (local_index >= residue.atom_count) {
            throw std::out_of_range("local atom index out of range");
        }
        return {molecule, residue.atom_begin + local_index};
    }

    AtomView name2atom(const std::string& name) const {
        const auto& residue = get();
        for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
            const AtomId atom_id = residue.atom_begin + local;
            if (molecule->atom(atom_id).name == name) {
                return {molecule, atom_id};
            }
        }
        throw std::out_of_range("atom name not found in residue: " + name);
    }
};

std::vector<ResidueView> residue_views(const std::shared_ptr<Molecule>& molecule) {
    std::vector<ResidueView> views;
    views.reserve(molecule->residue_count());
    for (ResidueId i = 0; i < molecule->residue_count(); ++i) {
        views.push_back({molecule, i});
    }
    return views;
}

std::vector<AtomView> residue_atom_views(const ResidueView& residue_view) {
    std::vector<AtomView> views;
    const auto& residue = residue_view.get();
    views.reserve(residue.atom_count);
    for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
        views.push_back({residue_view.molecule, residue.atom_begin + local});
    }
    return views;
}

std::shared_ptr<Molecule> load_pdb_object(py::object source) {
    return std::make_shared<Molecule>(load_pdb_text(read_python_input(source)));
}

std::shared_ptr<Molecule> load_mol2_object(py::object source) {
    return std::make_shared<Molecule>(load_mol2_text(read_python_input(source)));
}

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

std::shared_ptr<Molecule> get_template_molecule_object(const std::string& name) {
    return std::make_shared<Molecule>(get_template_molecule(name));
}

void set_box_padding_object(const std::shared_ptr<Molecule>& molecule, double padding, bool center) {
    molecule->set_box_padding(padding, center);
}

void add_solvent_box_object(const std::shared_ptr<Molecule>& molecule, const std::shared_ptr<Molecule>& solvent,
                            double distance, double tolerance, py::object n_solvent, std::uint64_t seed) {
    std::int64_t count = 0;
    if (!n_solvent.is_none()) {
        count = py::cast<std::int64_t>(n_solvent);
    }
    add_solvent_box(*molecule, *solvent, distance, tolerance, count, seed);
}

void add_ions_object(const std::shared_ptr<Molecule>& molecule,
                     const std::unordered_map<std::string, std::int64_t>& counts, std::uint64_t seed,
                     const std::string& solvent) {
    add_ions(*molecule, counts, seed, solvent);
}

void add_molecule_object(const std::shared_ptr<Molecule>& molecule, const std::shared_ptr<Molecule>& other) {
    molecule->add_molecule(*other);
}

std::unordered_map<std::string, std::string> save_sponge_input_object(const std::shared_ptr<Molecule>& molecule,
                                                                      const std::string& prefix,
                                                                      const std::string& dirname) {
    const auto paths = save_sponge_input(*molecule, prefix.empty() ? molecule->name : prefix, dirname);
    std::unordered_map<std::string, std::string> out;
    for (const auto& [key, value] : paths) {
        out[key] = value.string();
    }
    return out;
}

void save_pdb_object(const std::shared_ptr<Molecule>& molecule, const std::string& filename) {
    save_pdb(*molecule, filename);
}

void save_mol2_object(const std::shared_ptr<Molecule>& molecule, const std::string& filename) {
    save_mol2(*molecule, filename);
}

}  // namespace

PYBIND11_MODULE(_core, m) {
    py::class_<AtomView>(m, "Atom")
        .def_property_readonly("name", [](const AtomView& self) { return self.get().name; })
        .def_property_readonly("type", [](const AtomView& self) { return self.get().type; })
        .def_property_readonly("element", [](const AtomView& self) { return self.get().element; })
        .def_property("x", [](const AtomView& self) { return self.get().x; },
                      [](AtomView& self, double value) { self.molecule->atom(self.id).x = value; })
        .def_property("y", [](const AtomView& self) { return self.get().y; },
                      [](AtomView& self, double value) { self.molecule->atom(self.id).y = value; })
        .def_property("z", [](const AtomView& self) { return self.get().z; },
                      [](AtomView& self, double value) { self.molecule->atom(self.id).z = value; })
        .def_property_readonly("charge", [](const AtomView& self) { return self.get().charge; })
        .def_property_readonly("mass", [](const AtomView& self) { return self.get().mass; });

    py::class_<ResidueView>(m, "Residue")
        .def_property_readonly("name", [](const ResidueView& self) { return self.get().name; })
        .def_property_readonly("type_name", [](const ResidueView& self) { return self.get().type_name; })
        .def_property_readonly("atom_count", [](const ResidueView& self) { return self.get().atom_count; })
        .def_property_readonly("atoms", &residue_atom_views)
        .def("atom", &ResidueView::atom)
        .def("name2atom", &ResidueView::name2atom)
        .def("Name2Atom", &ResidueView::name2atom);

    py::class_<Molecule, std::shared_ptr<Molecule>>(m, "Molecule")
        .def(py::init<std::string>(), py::arg("name") = "MOL")
        .def_readwrite("name", &Molecule::name)
        .def_property_readonly("atom_count", &Molecule::atom_count)
        .def_property_readonly("residue_count", &Molecule::residue_count)
        .def_property_readonly("residues", &residue_views)
        .def_readwrite("box_length", &Molecule::box_length)
        .def_readwrite("box_angle", &Molecule::box_angle)
        .def("set_box_padding", &Molecule::set_box_padding, py::arg("padding") = 0.5, py::arg("center") = true)
        .def("Set_Box_Padding", &Molecule::set_box_padding, py::arg("padding") = 0.5, py::arg("center") = true)
        .def("add_molecule", &Molecule::add_molecule, py::arg("other"))
        .def("Add_Molecule", &Molecule::add_molecule, py::arg("other"))
        .def("validate", &Molecule::validate)
        .def("residue_counts", &Molecule::residue_counts);

    py::class_<ResidueType>(m, "ResidueType")
        .def(py::init<std::string>())
        .def_property_readonly("name", &ResidueType::name)
        .def_property_readonly("version", &ResidueType::version)
        .def_property_readonly("atom_count", &ResidueType::atom_count)
        .def_property_readonly("bond_count", &ResidueType::bond_count)
        .def("add_atom", &ResidueType::add_atom, py::arg("name"), py::arg("atom_type"), py::arg("x"),
             py::arg("y"), py::arg("z"), py::arg("charge") = 0.0, py::arg("mass") = 0.0)
        .def("Add_Atom", &ResidueType::add_atom, py::arg("name"), py::arg("atom_type"), py::arg("x"),
             py::arg("y"), py::arg("z"), py::arg("charge") = 0.0, py::arg("mass") = 0.0)
        .def("add_connectivity", &ResidueType::add_connectivity)
        .def("Add_Connectivity", &ResidueType::add_connectivity);

    py::class_<Assign, std::shared_ptr<Assign>>(m, "Assign")
        .def(py::init<std::string>(), py::arg("name") = "ASN")
        .def_readwrite("name", &Assign::name)
        .def_property_readonly("atom_count", &Assign::atom_count)
        .def_property_readonly("bond_count", &Assign::bond_count)
        .def_property_readonly("atoms", [](const Assign& self) { return self.elements; })
        .def_readonly("element_details", &Assign::element_details)
        .def_readonly("atom_types", &Assign::atom_types)
        .def("add_atom", &Assign::add_atom, py::arg("element"), py::arg("x"), py::arg("y"), py::arg("z"),
             py::arg("name") = "", py::arg("charge") = 0.0)
        .def("Add_Atom", &Assign::add_atom, py::arg("element"), py::arg("x"), py::arg("y"), py::arg("z"),
             py::arg("name") = "", py::arg("charge") = 0.0)
        .def("add_bond", &Assign::add_bond, py::arg("atom1"), py::arg("atom2"), py::arg("order") = 1)
        .def("determine_connectivity", &Assign::determine_connectivity, py::arg("simple_cutoff"))
        .def("determine_bond_order", &Assign::determine_bond_order, py::arg("check_formal_charge") = true,
             py::arg("total_charge") = py::none())
        .def("determine_atom_type", &Assign::determine_atom_type, py::arg("rule"))
        .def("to_residuetype", &Assign::to_residuetype, py::arg("name"))
        .def("to_molecule", [](const std::shared_ptr<Assign>& self, const std::string& name) {
            return std::make_shared<Molecule>(self->to_molecule(name));
        }, py::arg("name"));

    m.def("load_pdb", &load_pdb_object);
    m.def("load_mol2", &load_mol2_object);
    m.def("get_assignment_from_mol2", &get_assignment_from_mol2_object, py::arg("source"),
          py::arg("total_charge") = py::none());
    m.def("load_frcmod", [](const std::string& filename) {
        register_amber_frcmod_file(filename);
        return py::dict();
    });
    m.def("load_parmdat", [](const std::string& filename) {
        register_amber_parmdat_file(filename);
        return py::dict();
    });
    m.def("add_solvent_box", &add_solvent_box_object, py::arg("molecule"), py::arg("solvent"),
          py::arg("distance"), py::arg("tolerance") = 2.5, py::arg("n_solvent") = py::none(),
          py::arg("seed") = 0);
    m.def("add_ions", &add_ions_object, py::arg("molecule"), py::arg("counts"), py::arg("seed") = 0,
          py::arg("solvent") = "WAT");
    m.def("add_molecule", &add_molecule_object, py::arg("molecule"), py::arg("other"));
    m.def("set_box_padding", &set_box_padding_object, py::arg("molecule"), py::arg("padding") = 0.5,
          py::arg("center") = true);
    m.def("save_sponge_input", &save_sponge_input_object, py::arg("molecule"), py::arg("prefix") = "",
          py::arg("dirname") = ".");
    m.def("save_pdb", &save_pdb_object, py::arg("molecule"), py::arg("filename"));
    m.def("save_mol2", &save_mol2_object, py::arg("molecule"), py::arg("filename"));
    m.def("implemented_gaff_assign_types", &implemented_gaff_assign_types);

    m.def("register_ff14sb", &register_ff14sb);
    m.def("register_tip3p", &register_tip3p);
    m.def("register_amber_parmdat_file", [](const std::string& filename) { register_amber_parmdat_file(filename); });
    m.def("register_amber_frcmod_file", [](const std::string& filename) { register_amber_frcmod_file(filename); });
    m.def("register_residue_templates_from_mol2_file",
          [](const std::string& filename) { register_residue_templates_from_mol2_file(filename); });
    m.def("has_template", &has_template);
    m.def("template_atom_count", &template_atom_count);
    m.def("get_template_molecule", &get_template_molecule_object);
}
