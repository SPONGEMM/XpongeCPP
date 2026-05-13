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

std::vector<AtomView> molecule_atom_views(const std::shared_ptr<Molecule>& molecule) {
    std::vector<AtomView> views;
    views.reserve(molecule->atom_count());
    for (AtomId i = 0; i < molecule->atom_count(); ++i) {
        views.push_back({molecule, i});
    }
    return views;
}

std::vector<std::array<AtomId, 2>> molecule_explicit_bonds(const std::shared_ptr<Molecule>& molecule) {
    std::vector<std::array<AtomId, 2>> bonds;
    bonds.reserve(molecule->explicit_bonds.size());
    for (const auto& bond : molecule->explicit_bonds) {
        bonds.push_back({bond.atom1, bond.atom2});
    }
    return bonds;
}

std::vector<std::array<AtomId, 2>> molecule_residue_links(const std::shared_ptr<Molecule>& molecule) {
    std::vector<std::array<AtomId, 2>> links;
    links.reserve(molecule->residue_links.size());
    for (const auto& link : molecule->residue_links) {
        links.push_back({link.atom1, link.atom2});
    }
    return links;
}

std::shared_ptr<Molecule> load_pdb_object(py::object source, bool judge_histone, const std::string& position_need,
                                          bool ignore_hydrogen, bool ignore_unknown_name, bool ignore_seqres,
                                          bool ignore_conect, bool read_cryst1, py::object unterminal_residues) {
    PdbLoadOptions options;
    options.judge_histone = judge_histone;
    options.position_need = position_need.empty() ? 'A' : position_need[0];
    options.ignore_hydrogen = ignore_hydrogen;
    options.ignore_unknown_name = ignore_unknown_name;
    options.ignore_seqres = ignore_seqres;
    options.ignore_conect = ignore_conect;
    options.read_cryst1 = read_cryst1;
    if (!unterminal_residues.is_none()) {
        for (const auto item : unterminal_residues) {
            options.unterminal_residues.push_back(py::str(item));
        }
    }
    return std::make_shared<Molecule>(load_pdb_text(read_python_input(source), options));
}

std::shared_ptr<Molecule> load_mol2_object(py::object source) {
    return std::make_shared<Molecule>(load_mol2_text(read_python_input(source)));
}

std::shared_ptr<Molecule> load_gromacs_topology_object(const std::string& filename) {
    return std::make_shared<Molecule>(load_gromacs_topology_file(filename));
}

std::shared_ptr<Molecule> load_opls_itp_object(const std::string& filename) {
    return std::make_shared<Molecule>(load_opls_itp_file(filename));
}

std::shared_ptr<Molecule> load_charmm_topology_object(const std::string& filename) {
    return std::make_shared<Molecule>(load_charmm_topology_file(filename));
}

void load_sw_parameter_object(const std::string& filename, const std::shared_ptr<Molecule>& molecule) {
    load_sw_parameter_file(filename, *molecule);
}

void load_edip_parameter_object(const std::string& filename, const std::shared_ptr<Molecule>& molecule) {
    load_edip_parameter_file(filename, *molecule);
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

std::shared_ptr<Molecule> get_template_molecule_object(const std::string& name) {
    return std::make_shared<Molecule>(get_template_molecule(name));
}

void set_box_padding_object(const std::shared_ptr<Molecule>& molecule, double padding, bool center) {
    molecule->set_box_padding(padding, center);
}

std::array<double, 6> parse_solvent_distance(py::object distance) {
    if (PySequence_Check(distance.ptr()) != 0 && !PyUnicode_Check(distance.ptr()) && !PyBytes_Check(distance.ptr())) {
        const py::sequence values = py::reinterpret_borrow<py::sequence>(distance);
        if (values.size() == 3) {
            return {py::cast<double>(values[0]), py::cast<double>(values[1]), py::cast<double>(values[2]),
                    py::cast<double>(values[0]), py::cast<double>(values[1]), py::cast<double>(values[2])};
        }
        if (values.size() == 6) {
            return {py::cast<double>(values[0]), py::cast<double>(values[1]), py::cast<double>(values[2]),
                    py::cast<double>(values[3]), py::cast<double>(values[4]), py::cast<double>(values[5])};
        }
        throw py::type_error("the length of parameter distance should be 3 or 6");
    }
    const double scalar = py::cast<double>(distance);
    return {scalar, scalar, scalar, scalar, scalar, scalar};
}

void add_solvent_box_object(const std::shared_ptr<Molecule>& molecule, const std::shared_ptr<Molecule>& solvent,
                            py::object distance, double tolerance, py::object n_solvent, std::uint64_t seed) {
    std::int64_t count = 0;
    if (!n_solvent.is_none()) {
        count = py::cast<std::int64_t>(n_solvent);
    }
    add_solvent_box(*molecule, *solvent, parse_solvent_distance(distance), tolerance, count, seed);
}

void add_ions_object(const std::shared_ptr<Molecule>& molecule,
                     const std::unordered_map<std::string, std::int64_t>& counts, std::uint64_t seed,
                     const std::string& solvent) {
    add_ions(*molecule, counts, seed, solvent);
}

void add_molecule_object(const std::shared_ptr<Molecule>& molecule, const std::shared_ptr<Molecule>& other) {
    molecule->add_molecule(*other);
}

void replace_residues_object(const std::shared_ptr<Molecule>& molecule,
                             const std::unordered_map<ResidueId, std::shared_ptr<Molecule>>& replacements,
                             const std::vector<double>& residue_sort_keys, bool sort) {
    std::unordered_map<ResidueId, Molecule> copied_replacements;
    copied_replacements.reserve(replacements.size());
    for (const auto& [residue_id, replacement] : replacements) {
        copied_replacements.emplace(residue_id, *replacement);
    }
    molecule->replace_residues(copied_replacements, residue_sort_keys, sort);
}

void reorder_atoms_by_template_object(const std::shared_ptr<Molecule>& molecule,
                                      const std::shared_ptr<Molecule>& template_molecule) {
    molecule->reorder_atoms_by_template(*template_molecule);
}

std::shared_ptr<Molecule> molecule_binary_add(const std::shared_ptr<Molecule>& lhs,
                                              const std::shared_ptr<Molecule>& rhs, bool link) {
    auto out = std::make_shared<Molecule>(*lhs);
    out->add_molecule_linked(*rhs, link);
    return out;
}

std::shared_ptr<Molecule> molecule_repeat(const std::shared_ptr<Molecule>& molecule, int count) {
    if (count < 1) {
        throw std::invalid_argument("multiple should be not less than 1");
    }
    auto out = std::make_shared<Molecule>(*molecule);
    for (int i = 1; i < count; ++i) {
        out->add_molecule_linked(*molecule, true);
    }
    return out;
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

void save_pdb_object(const std::shared_ptr<Molecule>& molecule, const std::string& filename, bool write_cryst1) {
    if (write_cryst1) {
        save_pdb(*molecule, filename);
        return;
    }
    Molecule copy = *molecule;
    copy.has_box = false;
    save_pdb(copy, filename);
}

void save_mol2_object(const std::shared_ptr<Molecule>& molecule, const std::string& filename) {
    save_mol2(*molecule, filename);
}

py::tuple gro_data_to_python(const GroData& data) {
    std::vector<std::vector<double>> coordinates;
    coordinates.reserve(data.coordinates.size());
    for (const auto& coordinate : data.coordinates) {
        coordinates.push_back({coordinate[0], coordinate[1], coordinate[2]});
    }
    std::vector<double> box{data.box[0], data.box[1], data.box[2], data.box[3], data.box[4], data.box[5]};
    return py::make_tuple(coordinates, box);
}

py::tuple coordinate_data_to_python(const CoordinateData& data) {
    std::vector<std::vector<double>> coordinates;
    coordinates.reserve(data.coordinates.size());
    for (const auto& coordinate : data.coordinates) {
        coordinates.push_back({coordinate[0], coordinate[1], coordinate[2]});
    }
    std::vector<double> box;
    if (data.has_box) {
        box = {data.box[0], data.box[1], data.box[2], data.box[3], data.box[4], data.box[5]};
    }
    return py::make_tuple(coordinates, box);
}

py::tuple load_coordinate_object(py::object source, py::object molecule) {
    if (molecule.is_none()) {
        return coordinate_data_to_python(load_coordinate_text(read_python_input(source)));
    }
    auto mol = py::cast<std::shared_ptr<Molecule>>(molecule);
    return coordinate_data_to_python(load_coordinate_text(read_python_input(source), *mol));
}

py::tuple load_rst7_object(py::object source, py::object molecule) {
    if (molecule.is_none()) {
        return coordinate_data_to_python(load_rst7_text(read_python_input(source)));
    }
    auto mol = py::cast<std::shared_ptr<Molecule>>(molecule);
    return coordinate_data_to_python(load_rst7_text(read_python_input(source), *mol));
}

py::tuple load_gro_object(py::object source, py::object molecule, bool read_box_angle) {
    if (molecule.is_none()) {
        return gro_data_to_python(load_gro_text(read_python_input(source), read_box_angle));
    }
    auto mol = py::cast<std::shared_ptr<Molecule>>(molecule);
    return gro_data_to_python(load_gro_text(read_python_input(source), *mol, read_box_angle));
}

void save_gro_object(const std::shared_ptr<Molecule>& molecule, const std::string& filename) {
    save_gro(*molecule, filename);
}

py::tuple load_molpsf_object(py::object source, py::object split_by) {
    std::string split = "connectivity";
    if (!split_by.is_none()) {
        split = py::cast<std::string>(split_by);
    } else {
        split.clear();
    }
    auto data = load_molpsf_text(read_python_input(source), split);
    py::dict molecules;
    for (auto& [name, molecule] : data.molecules) {
        molecules[py::str(name)] = std::make_shared<Molecule>(std::move(molecule));
    }
    return py::make_tuple(std::make_shared<Molecule>(std::move(data.molecule)), molecules);
}

py::tuple merge_dual_topology_object(const std::shared_ptr<Molecule>& molecule, ResidueId residue_index,
                                     const std::shared_ptr<Molecule>& residue_b_molecule,
                                     const std::unordered_map<std::uint32_t, std::uint32_t>& match_b_to_a) {
    auto [from_a, to_b] = merge_dual_topology(*molecule, residue_index, *residue_b_molecule, match_b_to_a);
    return py::make_tuple(std::make_shared<Molecule>(std::move(from_a)),
                          std::make_shared<Molecule>(std::move(to_b)),
                          match_b_to_a);
}

std::shared_ptr<Molecule> merge_force_field_object(
    const std::shared_ptr<Molecule>& molecule_a, const std::shared_ptr<Molecule>& molecule_b, double default_lambda,
    const std::unordered_map<std::string, double>& specific_lambda) {
    return std::make_shared<Molecule>(merge_force_field(*molecule_a, *molecule_b, default_lambda, specific_lambda));
}

}  // namespace

PYBIND11_MODULE(_core, m) {
    py::class_<AtomView>(m, "Atom")
        .def_property_readonly("index", [](const AtomView& self) { return self.id; })
        .def_property_readonly("name", [](const AtomView& self) { return self.get().name; })
        .def_property("type", [](const AtomView& self) { return self.get().type; },
                      [](AtomView& self, const std::string& value) { self.molecule->atom(self.id).type = value; })
        .def_property_readonly("element", [](const AtomView& self) { return self.get().element; })
        .def_property("x", [](const AtomView& self) { return self.get().x; },
                      [](AtomView& self, double value) { self.molecule->atom(self.id).x = value; })
        .def_property("y", [](const AtomView& self) { return self.get().y; },
                      [](AtomView& self, double value) { self.molecule->atom(self.id).y = value; })
        .def_property("z", [](const AtomView& self) { return self.get().z; },
                      [](AtomView& self, double value) { self.molecule->atom(self.id).z = value; })
        .def_property("charge", [](const AtomView& self) { return self.get().charge; },
                      [](AtomView& self, double value) { self.molecule->atom(self.id).charge = value; })
        .def_property("mass", [](const AtomView& self) { return self.get().mass; },
                      [](AtomView& self, double value) { self.molecule->atom(self.id).mass = value; })
        .def_property("lj_type_b", [](const AtomView& self) { return self.get().lj_type_b; },
                      [](AtomView& self, const std::string& value) { self.molecule->atom(self.id).lj_type_b = value; })
        .def_property("sw_type", [](const AtomView& self) { return self.get().sw_type; },
                      [](AtomView& self, const std::string& value) { self.molecule->atom(self.id).sw_type = value; })
        .def_property("edip_type", [](const AtomView& self) { return self.get().edip_type; },
                      [](AtomView& self, const std::string& value) { self.molecule->atom(self.id).edip_type = value; })
        .def_property("gb_radius", [](const AtomView& self) { return self.get().gb_radius; },
                      [](AtomView& self, double value) { self.molecule->atom(self.id).gb_radius = value; })
        .def_property("gb_scaler", [](const AtomView& self) { return self.get().gb_scaler; },
                      [](AtomView& self, double value) { self.molecule->atom(self.id).gb_scaler = value; })
        .def_property("subsys", [](const AtomView& self) { return self.get().subsys; },
                      [](AtomView& self, int value) { self.molecule->atom(self.id).subsys = value; })
        .def_property("bad_coordinate", [](const AtomView& self) { return self.get().bad_coordinate; },
                      [](AtomView& self, bool value) { self.molecule->atom(self.id).bad_coordinate = value; })
        .def_property("zero_lj_atom", [](const AtomView& self) { return self.get().zero_lj_atom; },
                      [](AtomView& self, bool value) { self.molecule->atom(self.id).zero_lj_atom = value; })
        .def_property_readonly("serial", [](const AtomView& self) { return self.get().serial; })
        .def_property_readonly("altloc", [](const AtomView& self) { return std::string(1, self.get().altloc); })
        .def_property_readonly("occupancy", [](const AtomView& self) { return self.get().occupancy; })
        .def_property_readonly("temp_factor", [](const AtomView& self) { return self.get().temp_factor; })
        .def_property_readonly("record_name", [](const AtomView& self) { return self.get().record_name; });

    py::class_<ResidueView>(m, "Residue")
        .def_property_readonly("index", [](const ResidueView& self) { return self.id; })
        .def_property_readonly("name", [](const ResidueView& self) { return self.get().name; })
        .def_property_readonly("type_name", [](const ResidueView& self) { return self.get().type_name; })
        .def_property_readonly("chain_id", [](const ResidueView& self) { return std::string(1, self.get().chain_id); })
        .def_property_readonly("effective_chain_id",
                               [](const ResidueView& self) { return std::string(1, self.get().effective_chain_id); })
        .def_property_readonly("segment_id", [](const ResidueView& self) { return self.get().segment_id; })
        .def_property_readonly("pdb_resseq", [](const ResidueView& self) { return self.get().pdb_resseq; })
        .def_property_readonly("insertion_code",
                               [](const ResidueView& self) { return std::string(1, self.get().insertion_code); })
        .def_property_readonly("is_hetero", [](const ResidueView& self) { return self.get().is_hetero; })
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
        .def_property_readonly("atoms", &molecule_atom_views)
        .def_property_readonly("residues", &residue_views)
        .def_property_readonly("explicit_bonds", &molecule_explicit_bonds)
        .def_property_readonly("residue_links", &molecule_residue_links)
        .def_readwrite("box_length", &Molecule::box_length)
        .def_readwrite("box_angle", &Molecule::box_angle)
        .def("set_box_padding", &Molecule::set_box_padding, py::arg("padding") = 0.5, py::arg("center") = true)
        .def("Set_Box_Padding", &Molecule::set_box_padding, py::arg("padding") = 0.5, py::arg("center") = true)
        .def("add_molecule", &Molecule::add_molecule, py::arg("other"))
        .def("Add_Molecule", &Molecule::add_molecule, py::arg("other"))
        .def("add_residue_link", &Molecule::add_residue_link, py::arg("atom1"), py::arg("atom2"))
        .def("Add_Residue_Link", &Molecule::add_residue_link, py::arg("atom1"), py::arg("atom2"))
        .def("copy", [](const std::shared_ptr<Molecule>& self) { return std::make_shared<Molecule>(*self); })
        .def("deepcopy", [](const std::shared_ptr<Molecule>& self) { return std::make_shared<Molecule>(*self); })
        .def("__add__", [](const std::shared_ptr<Molecule>& self, const std::shared_ptr<Molecule>& other) {
            return molecule_binary_add(self, other, true);
        }, py::is_operator())
        .def("__or__", [](const std::shared_ptr<Molecule>& self, const std::shared_ptr<Molecule>& other) {
            return molecule_binary_add(self, other, false);
        }, py::is_operator())
        .def("__mul__", [](const std::shared_ptr<Molecule>& self, int count) { return molecule_repeat(self, count); },
             py::is_operator())
        .def("__rmul__", [](const std::shared_ptr<Molecule>& self, int count) { return molecule_repeat(self, count); },
             py::is_operator())
        .def("__iadd__", [](const std::shared_ptr<Molecule>& self, const std::shared_ptr<Molecule>& other) {
            self->add_molecule_linked(*other, true);
            return self;
        }, py::is_operator())
        .def("__ior__", [](const std::shared_ptr<Molecule>& self, const std::shared_ptr<Molecule>& other) {
            self->add_molecule_linked(*other, false);
            return self;
        }, py::is_operator())
        .def("__imul__", [](const std::shared_ptr<Molecule>& self, int count) {
            if (count < 1) {
                throw std::invalid_argument("multiple should be not less than 1");
            }
            const Molecule base = *self;
            for (int i = 1; i < count; ++i) {
                self->add_molecule_linked(base, true);
            }
            return self;
        }, py::is_operator())
        .def("add_virtual_atom2", &Molecule::add_virtual_atom2, py::arg("virtual_atom"), py::arg("atom0"),
             py::arg("atom1"), py::arg("atom2"), py::arg("k1"), py::arg("k2"))
        .def("Add_Virtual_Atom2", &Molecule::add_virtual_atom2, py::arg("virtual_atom"), py::arg("atom0"),
             py::arg("atom1"), py::arg("atom2"), py::arg("k1"), py::arg("k2"))
        .def("add_improper_dihedral", &Molecule::add_improper_dihedral, py::arg("atom0"), py::arg("atom1"),
             py::arg("atom2"), py::arg("atom3"), py::arg("k"), py::arg("phi0"))
        .def("Add_Improper_Dihedral", &Molecule::add_improper_dihedral, py::arg("atom0"), py::arg("atom1"),
             py::arg("atom2"), py::arg("atom3"), py::arg("k"), py::arg("phi0"))
        .def("add_cmap_type", &Molecule::add_cmap_type, py::arg("resolution"), py::arg("parameters"))
        .def("Add_CMap_Type", &Molecule::add_cmap_type, py::arg("resolution"), py::arg("parameters"))
        .def("add_cmap", &Molecule::add_cmap, py::arg("atom0"), py::arg("atom1"), py::arg("atom2"),
             py::arg("atom3"), py::arg("atom4"), py::arg("type"))
        .def("Add_CMap", &Molecule::add_cmap, py::arg("atom0"), py::arg("atom1"), py::arg("atom2"),
             py::arg("atom3"), py::arg("atom4"), py::arg("type"))
        .def("add_nb14_extra", &Molecule::add_nb14_extra, py::arg("atom1"), py::arg("atom2"), py::arg("a"),
             py::arg("b"), py::arg("kee"))
        .def("Add_NB14_Extra", &Molecule::add_nb14_extra, py::arg("atom1"), py::arg("atom2"), py::arg("a"),
             py::arg("b"), py::arg("kee"))
        .def("add_urey_bradley", &Molecule::add_urey_bradley, py::arg("atom0"), py::arg("atom1"),
             py::arg("atom2"), py::arg("k"), py::arg("b"), py::arg("kUB"), py::arg("r13"))
        .def("Add_Urey_Bradley", &Molecule::add_urey_bradley, py::arg("atom0"), py::arg("atom1"),
             py::arg("atom2"), py::arg("k"), py::arg("b"), py::arg("kUB"), py::arg("r13"))
        .def("add_ryckaert_bellemans", &Molecule::add_ryckaert_bellemans, py::arg("atom0"),
             py::arg("atom1"), py::arg("atom2"), py::arg("atom3"), py::arg("c0"), py::arg("c1"),
             py::arg("c2"), py::arg("c3"), py::arg("c4"), py::arg("c5"))
        .def("Add_Ryckaert_Bellemans", &Molecule::add_ryckaert_bellemans, py::arg("atom0"),
             py::arg("atom1"), py::arg("atom2"), py::arg("atom3"), py::arg("c0"), py::arg("c1"),
             py::arg("c2"), py::arg("c3"), py::arg("c4"), py::arg("c5"))
        .def("add_bond_soft", &Molecule::add_bond_soft, py::arg("atom1"), py::arg("atom2"),
             py::arg("k"), py::arg("b"), py::arg("from_AorB"))
        .def("Add_Bond_Soft", &Molecule::add_bond_soft, py::arg("atom1"), py::arg("atom2"),
             py::arg("k"), py::arg("b"), py::arg("from_AorB"))
        .def("add_listed_force_definition", &Molecule::add_listed_force_definition, py::arg("definition"))
        .def("Add_Listed_Force_Definition", &Molecule::add_listed_force_definition, py::arg("definition"))
        .def("set_gb_radius", &Molecule::set_gb_radius, py::arg("radius_set") = "modified_bondi_radii")
        .def("Set_GB_Radius", &Molecule::set_gb_radius, py::arg("radius_set") = "modified_bondi_radii")
        .def("enable_min_bonded_parameters", &Molecule::enable_min_bonded_parameters, py::arg("enabled") = true)
        .def("Enable_Min_Bonded_Parameters", &Molecule::enable_min_bonded_parameters, py::arg("enabled") = true)
        .def("enable_subsys_division", &Molecule::enable_subsys_division, py::arg("enabled") = true)
        .def("Enable_Subsys_Division", &Molecule::enable_subsys_division, py::arg("enabled") = true)
        .def("enable_lj_soft_core", &Molecule::enable_lj_soft_core, py::arg("enabled") = true)
        .def("Enable_LJ_Soft_Core", &Molecule::enable_lj_soft_core, py::arg("enabled") = true)
        .def("add_sw_type", &Molecule::add_sw_type, py::arg("name"), py::arg("A"), py::arg("B"),
             py::arg("epsilon"), py::arg("p"), py::arg("q"), py::arg("a"), py::arg("gamma"),
             py::arg("sigma"), py::arg("lambda_"), py::arg("b"))
        .def("Add_SW_Type", &Molecule::add_sw_type, py::arg("name"), py::arg("A"), py::arg("B"),
             py::arg("epsilon"), py::arg("p"), py::arg("q"), py::arg("a"), py::arg("gamma"),
             py::arg("sigma"), py::arg("lambda_"), py::arg("b"))
        .def("add_edip_type", &Molecule::add_edip_type, py::arg("name"), py::arg("A"), py::arg("B"),
             py::arg("a"), py::arg("c"), py::arg("alpha"), py::arg("beta"), py::arg("eta"),
             py::arg("gamma"), py::arg("lambda_"), py::arg("mu"), py::arg("rho"), py::arg("sigma"),
             py::arg("Q0"), py::arg("u1"), py::arg("u2"), py::arg("u3"), py::arg("u4"))
        .def("Add_EDIP_Type", &Molecule::add_edip_type, py::arg("name"), py::arg("A"), py::arg("B"),
             py::arg("a"), py::arg("c"), py::arg("alpha"), py::arg("beta"), py::arg("eta"),
             py::arg("gamma"), py::arg("lambda_"), py::arg("mu"), py::arg("rho"), py::arg("sigma"),
             py::arg("Q0"), py::arg("u1"), py::arg("u2"), py::arg("u3"), py::arg("u4"))
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
        .def_readonly("names", &Assign::names)
        .def_readonly("coordinates", &Assign::coordinates)
        .def_readonly("charges", &Assign::charges)
        .def_readonly("formal_charges", &Assign::formal_charges)
        .def_readonly("bonds", &Assign::bonds)
        .def_readonly("atom_types", &Assign::atom_types)
        .def("add_atom", &Assign::add_atom, py::arg("element"), py::arg("x"), py::arg("y"), py::arg("z"),
             py::arg("name") = "", py::arg("charge") = 0.0)
        .def("Add_Atom", &Assign::add_atom, py::arg("element"), py::arg("x"), py::arg("y"), py::arg("z"),
             py::arg("name") = "", py::arg("charge") = 0.0)
        .def("add_bond", &Assign::add_bond, py::arg("atom1"), py::arg("atom2"), py::arg("order") = 1)
        .def("delete_atom", &Assign::delete_atom, py::arg("atom"))
        .def("set_charge", &Assign::set_charge, py::arg("atom"), py::arg("charge"))
        .def("set_charges", &Assign::set_charges, py::arg("charges"))
        .def("set_formal_charge", &Assign::set_formal_charge, py::arg("atom"), py::arg("charge"))
        .def("set_coordinate", &Assign::set_coordinate, py::arg("atom"), py::arg("x"), py::arg("y"), py::arg("z"))
        .def("set_atom_type", &Assign::set_atom_type, py::arg("atom"), py::arg("atom_type"))
        .def("has_atom_marker", &Assign::has_atom_marker, py::arg("atom"), py::arg("marker"))
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
             py::arg("check_formal_charge"),
             py::arg("total_charge"),
             py::arg("max_step"),
             py::arg("max_stat"),
             py::arg("penalty_scores"),
             py::arg("extra_criteria") = py::none())
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

    m.def("load_pdb", &load_pdb_object, py::arg("source"), py::arg("judge_histone") = true,
          py::arg("position_need") = "A", py::arg("ignore_hydrogen") = false,
          py::arg("ignore_unknown_name") = false, py::arg("ignore_seqres") = true,
          py::arg("ignore_conect") = true, py::arg("read_cryst1") = true,
          py::arg("unterminal_residues") = py::none());
    m.def("load_mol2", &load_mol2_object);
    m.def("load_molpsf", &load_molpsf_object, py::arg("source"), py::arg("split_by") = "connectivity");
    m.def("load_gromacs_topology_file", &load_gromacs_topology_object);
    m.def("load_opls_itp_file", &load_opls_itp_object);
    m.def("load_charmm_parameter_file", [](const std::string& filename) { load_charmm_parameter_file(filename); });
    m.def("load_charmm_topology_file", &load_charmm_topology_object);
    m.def("load_sw_parameter_file", &load_sw_parameter_object, py::arg("filename"), py::arg("molecule"));
    m.def("load_edip_parameter_file", &load_edip_parameter_object, py::arg("filename"), py::arg("molecule"));
    m.def("get_assignment_from_mol2", &get_assignment_from_mol2_object, py::arg("source"),
          py::arg("total_charge") = py::none());
    m.def("get_assignment_from_xyz", &get_assignment_from_xyz_object, py::arg("source"));
    m.def("get_assignment_from_pdb", &get_assignment_from_pdb_object, py::arg("source"));
    m.def("get_assignment_from_residuetype",
          [](const ResidueType& residue_type) {
              return std::make_shared<Assign>(get_assignment_from_residuetype(residue_type));
          }, py::arg("residue_type"));
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
    m.def("replace_residues", &replace_residues_object, py::arg("molecule"), py::arg("replacements"),
          py::arg("residue_sort_keys") = std::vector<double>{}, py::arg("sort") = true);
    m.def("reorder_atoms_by_template", &reorder_atoms_by_template_object, py::arg("molecule"),
          py::arg("template_molecule"));
    m.def("set_box_padding", &set_box_padding_object, py::arg("molecule"), py::arg("padding") = 0.5,
          py::arg("center") = true);
    m.def("save_sponge_input", &save_sponge_input_object, py::arg("molecule"), py::arg("prefix") = "",
          py::arg("dirname") = ".");
    m.def("save_pdb", &save_pdb_object, py::arg("molecule"), py::arg("filename"),
          py::arg("write_cryst1") = true);
    m.def("save_mol2", &save_mol2_object, py::arg("molecule"), py::arg("filename"));
    m.def("load_coordinate", &load_coordinate_object, py::arg("source"), py::arg("molecule") = py::none());
    m.def("load_rst7", &load_rst7_object, py::arg("source"), py::arg("molecule") = py::none());
    m.def("load_gro", &load_gro_object, py::arg("source"), py::arg("molecule") = py::none(),
          py::arg("read_box_angle") = true);
    m.def("save_gro", &save_gro_object, py::arg("molecule"), py::arg("filename"));
    m.def("implemented_gaff_assign_types", &implemented_gaff_assign_types);
    m.def("merge_dual_topology", &merge_dual_topology_object, py::arg("molecule"), py::arg("residue_index"),
          py::arg("residue_b_molecule"), py::arg("match_b_to_a"));
    m.def("merge_force_field", &merge_force_field_object, py::arg("molecule_a"), py::arg("molecule_b"),
          py::arg("default_lambda"), py::arg("specific_lambda") = std::unordered_map<std::string, double>{});

    m.def("register_ff14sb", &register_ff14sb);
    m.def("register_tip3p", &register_tip3p);
    m.def("register_amber_parmdat_file", [](const std::string& filename) { register_amber_parmdat_file(filename); });
    m.def("register_amber_frcmod_file", [](const std::string& filename) { register_amber_frcmod_file(filename); });
    m.def("register_amber_lj_parameter", &register_amber_lj_parameter, py::arg("atom_type"), py::arg("lj_type"),
          py::arg("epsilon"), py::arg("rmin"));
    m.def("register_amber_bond_parameter", &register_amber_bond_parameter, py::arg("atom_type1"),
          py::arg("atom_type2"), py::arg("k"), py::arg("length"));
    m.def("register_amber_cmap_parameter", &register_amber_cmap_parameter, py::arg("key"),
          py::arg("resolution"), py::arg("parameters"));
    m.def("register_residue_templates_from_mol2_text", &register_residue_templates_from_mol2_text,
          py::arg("text"));
    m.def("register_residue_templates_from_mol2_file",
          [](const std::string& filename) { register_residue_templates_from_mol2_file(filename); });
    m.def("register_template_molecule_from_mol2_file",
          [](const std::string& filename) { register_template_molecule_from_mol2_file(filename); });
    m.def("register_template_virtual_atom2", &register_template_virtual_atom2, py::arg("template_name"),
          py::arg("virtual_atom"), py::arg("atom0"), py::arg("atom1"), py::arg("atom2"), py::arg("k1"),
          py::arg("k2"));
    m.def("configure_residue_template_head", &configure_residue_template_head, py::arg("template_name"),
          py::arg("atom"), py::arg("length") = 1.5, py::arg("next") = "");
    m.def("configure_residue_template_tail", &configure_residue_template_tail, py::arg("template_name"),
          py::arg("atom"), py::arg("length") = 1.5, py::arg("next") = "");
    m.def("configure_residue_template_connect_atom", &configure_residue_template_connect_atom,
          py::arg("template_name"), py::arg("key"), py::arg("atom"));
    m.def("register_residue_template_alias", &register_residue_template_alias, py::arg("alias_name"),
          py::arg("template_name"));
    m.def("register_pdb_residue_name_mapping", &register_pdb_residue_name_mapping, py::arg("place"),
          py::arg("pdb_name"), py::arg("real_name"));
    m.def("register_pdb_residue_alias_mapping", &register_pdb_residue_alias_mapping, py::arg("pdb_name"),
          py::arg("real_name"));
    m.def("register_his_mapping", &register_his_mapping, py::arg("residue_name"), py::arg("hid"),
          py::arg("hie"), py::arg("hip"));
    m.def("has_template", &has_template);
    m.def("template_atom_count", &template_atom_count);
    m.def("get_template_molecule", &get_template_molecule_object);
    m.def("molecule_from_residuetype", [](const ResidueType& residue_type) {
        auto molecule = std::make_shared<Molecule>(residue_type.name());
        molecule->append_residue_from_type(residue_type, 0.0, 0.0, 0.0);
        return molecule;
    }, py::arg("residue_type"));
}
