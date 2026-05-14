#include "bindings_internal.hpp"

namespace xpongecpp {
namespace {

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

}  // namespace

void bind_forcefield_module(py::module_& m) {
    m.def("load_gromacs_topology_file", &load_gromacs_topology_object);
    m.def("load_opls_itp_file", &load_opls_itp_object);
    m.def("load_charmm_parameter_file", [](const std::string& filename) { load_charmm_parameter_file(filename); });
    m.def("load_charmm_topology_file", &load_charmm_topology_object);
    m.def("load_sw_parameter_file", &load_sw_parameter_object, py::arg("filename"), py::arg("molecule"));
    m.def("load_edip_parameter_file", &load_edip_parameter_object, py::arg("filename"), py::arg("molecule"));
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
    m.def("registered_template_names", &registered_template_names);
}

}  // namespace xpongecpp
