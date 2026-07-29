#include "bindings_internal.hpp"
#include "amber_internal.hpp"
#include "nonamber_internal.hpp"
#include "pdb_internal.hpp"

#include <mutex>

namespace xpongecpp {
namespace {

struct ForceFieldRegistrySnapshot {
    std::unordered_map<std::string, ResidueType> templates;
    std::unordered_map<std::string, Molecule> molecule_templates;
    std::vector<BondParameter> bond_parameters;
    std::vector<AngleParameter> angle_parameters;
    std::vector<ProperParameter> proper_parameters;
    std::vector<ImproperParameter> improper_parameters;
    std::vector<NB14Parameter> nb14_parameters;
    std::unordered_map<std::string, AmberCMapParameter> amber_cmap_parameters;
    std::unordered_map<std::string, std::string> lj_type_by_atom_type;
    std::unordered_map<std::string, std::pair<double, double>> lj_parameters;
    std::unordered_map<std::string, double> mass_by_atom_type;
    LJCombiningRule lj_combining_rule{LJCombiningRule::LorentzBerthelot};
    std::unordered_map<std::string, AtomTypeInfo> nonamber_atom_types;
    std::vector<CharmmUreyParameter> charmm_urey_parameters;
    std::vector<NBfixParameter> nbfix_parameters;
    std::unordered_map<std::string, std::string> pdb_head_map;
    std::unordered_map<std::string, std::string> pdb_tail_map;
    std::unordered_map<std::string, std::string> pdb_save_map;
    std::unordered_map<std::string, std::string> pdb_alias_map;
    std::unordered_map<std::string, HisNames> his_map;
};

ForceFieldRegistrySnapshot capture_forcefield_registries() {
    std::shared_lock lock(registry_mutex());
    return {
        templates(),
        molecule_templates(),
        bond_parameters(),
        angle_parameters(),
        proper_parameters(),
        improper_parameters(),
        nb14_parameters(),
        amber_cmap_parameters(),
        lj_type_by_atom_type(),
        lj_parameters(),
        mass_by_atom_type(),
        lj_combining_rule(),
        atom_type_registry(),
        charmm_urey_registry(),
        nbfix_registry(),
        pdb_head_map(),
        pdb_tail_map(),
        pdb_save_map(),
        pdb_alias_map(),
        his_map(),
    };
}

void restore_forcefield_registries(const ForceFieldRegistrySnapshot& snapshot) {
    std::unique_lock lock(registry_mutex());
    templates() = snapshot.templates;
    molecule_templates() = snapshot.molecule_templates;
    bond_parameters() = snapshot.bond_parameters;
    angle_parameters() = snapshot.angle_parameters;
    proper_parameters() = snapshot.proper_parameters;
    improper_parameters() = snapshot.improper_parameters;
    nb14_parameters() = snapshot.nb14_parameters;
    amber_cmap_parameters() = snapshot.amber_cmap_parameters;
    lj_type_by_atom_type() = snapshot.lj_type_by_atom_type;
    lj_parameters() = snapshot.lj_parameters;
    mass_by_atom_type() = snapshot.mass_by_atom_type;
    set_lj_combining_rule(snapshot.lj_combining_rule);
    atom_type_registry() = snapshot.nonamber_atom_types;
    charmm_urey_registry() = snapshot.charmm_urey_parameters;
    nbfix_registry() = snapshot.nbfix_parameters;
    pdb_head_map() = snapshot.pdb_head_map;
    pdb_tail_map() = snapshot.pdb_tail_map;
    pdb_save_map() = snapshot.pdb_save_map;
    pdb_alias_map() = snapshot.pdb_alias_map;
    his_map() = snapshot.his_map;
}

LJCombiningRule lj_rule_from_name(std::string name) {
    std::transform(name.begin(), name.end(), name.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    if (name == "lorentz_berthelot" || name == "lorentz-berthelot" || name == "lb") {
        return LJCombiningRule::LorentzBerthelot;
    }
    if (name == "good_hope" || name == "good-hope" || name == "goodhope") {
        return LJCombiningRule::GoodHope;
    }
    throw py::value_error("unsupported LJ combining rule: " + name);
}

std::string lj_rule_name(LJCombiningRule rule) {
    return rule == LJCombiningRule::GoodHope ? "good_hope" : "lorentz_berthelot";
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
    m.def("_snapshot_forcefield_registries", []() {
        auto* snapshot = new ForceFieldRegistrySnapshot(capture_forcefield_registries());
        return py::capsule(snapshot, "XpongeCPP.ForceFieldRegistrySnapshot", [](PyObject* capsule) {
            delete reinterpret_cast<ForceFieldRegistrySnapshot*>(
                PyCapsule_GetPointer(capsule, "XpongeCPP.ForceFieldRegistrySnapshot"));
        });
    });
    m.def("_restore_forcefield_registries", [](const py::capsule& capsule) {
        const auto* snapshot = reinterpret_cast<const ForceFieldRegistrySnapshot*>(
            capsule.get_pointer());
        if (snapshot == nullptr) {
            throw py::value_error("invalid force-field registry snapshot");
        }
        restore_forcefield_registries(*snapshot);
    }, py::arg("snapshot"));
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
    m.def("register_amber_angle_parameter", &register_amber_angle_parameter, py::arg("atom_types"),
          py::arg("k"), py::arg("theta"));
    m.def("register_amber_proper_dihedral_parameter", &register_amber_proper_dihedral_parameter,
          py::arg("atom_types"), py::arg("periodicity"), py::arg("k"), py::arg("phase"),
          py::arg("reset") = false);
    m.def("register_amber_improper_dihedral_parameter", &register_amber_improper_dihedral_parameter,
          py::arg("atom_types"), py::arg("periodicity"), py::arg("k"), py::arg("phase"));
    m.def("register_amber_nb14_scale", &register_amber_nb14_scale, py::arg("atom_type1"), py::arg("atom_type4"),
          py::arg("k_lj"), py::arg("k_ee"));
    m.def("register_amber_cmap_parameter", &register_amber_cmap_parameter, py::arg("key"),
          py::arg("resolution"), py::arg("parameters"));
    m.def("clear_amber_dihedral_parameters", &clear_amber_dihedral_parameters);
    m.def("clear_amber_improper_parameters", &clear_amber_improper_parameters);
    m.def("set_lj_combining_rule", [](const std::string& rule) {
        set_lj_combining_rule(lj_rule_from_name(rule));
    }, py::arg("rule"));
    m.def("get_lj_combining_rule", []() { return lj_rule_name(lj_combining_rule()); });
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
