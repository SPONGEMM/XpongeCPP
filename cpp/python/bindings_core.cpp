#include "bindings_internal.hpp"

#include <algorithm>
#include <initializer_list>

namespace xpongecpp {
namespace {

std::shared_ptr<Molecule> get_template_molecule_object(const std::string& name) {
    return std::make_shared<Molecule>(get_template_molecule(name));
}

void set_box_padding_object(const std::shared_ptr<Molecule>& molecule, double padding, bool center) {
    molecule->set_box_padding(padding, center);
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

void require_empty_native_protocol(const py::object& protocol) {
    if (!protocol.is_none()) {
        throw std::invalid_argument("nonempty native protocols are not supported yet");
    }
}

void require_native_bundle_compatibility(const std::shared_ptr<Molecule>& molecule) {
    const py::object workflows = py::module_::import("XpongeCPP._compat.workflows");
    if (py::cast<bool>(workflows.attr("min_bonded_parameters_enabled")())) {
        throw std::invalid_argument(
            "bundle export does not support minimum-bonded parameters "
            "(fake_mass, fake_LJ, fake_charge)");
    }

    const py::object molecule_type = py::module_::import("XpongeCPP._core").attr("Molecule");
    if (!py::hasattr(molecule_type, "_save_functions")) {
        return;
    }
    const py::dict serializers =
        py::reinterpret_borrow<py::dict>(molecule_type.attr("_save_functions"));
    if (py::len(serializers) == 0) {
        return;
    }
    py::list active_keys;
    for (const auto item : serializers) {
        const py::object payload =
            py::reinterpret_borrow<py::object>(item.second)(py::cast(molecule));
        const int is_active = PyObject_IsTrue(payload.ptr());
        if (is_active < 0) {
            throw py::error_already_set();
        }
        if (is_active != 0) {
            active_keys.append(item.first);
        }
    }
    if (active_keys.empty()) {
        return;
    }
    throw std::invalid_argument(
        "bundle export does not support active compatibility serializers " +
        py::str(active_keys).cast<std::string>() +
        "; native bundle conversion is required");
}

std::string normalize_bundle_prefix(const py::object& prefix, const std::string& fallback) {
    if (prefix.is_none()) {
        return fallback;
    }
    const std::string normalized = py::str(prefix).cast<std::string>();
    return normalized.empty() ? fallback : normalized;
}

std::shared_ptr<Molecule> save_sponge_input_bundle_object(
    const std::shared_ptr<Molecule>& molecule, py::object prefix, py::object dirname, py::object protocol) {
    require_empty_native_protocol(protocol);
    require_native_bundle_compatibility(molecule);
    save_sponge_input_bundle(*molecule, normalize_bundle_prefix(prefix, molecule->name),
                             py::str(dirname).cast<std::string>());
    return molecule;
}

std::shared_ptr<Molecule> save_residuetype_bundle_object(
    const ResidueType& residue_type, py::object prefix, py::object dirname, py::object protocol) {
    require_empty_native_protocol(protocol);
    auto molecule = std::make_shared<Molecule>(residue_type.name());
    molecule->append_residue_from_type(residue_type, 0.0, 0.0, 0.0);
    require_native_bundle_compatibility(molecule);
    save_sponge_input_bundle(*molecule, normalize_bundle_prefix(prefix, molecule->name),
                             py::str(dirname).cast<std::string>());
    return molecule;
}

std::shared_ptr<Molecule> molecule_from_residue_view(const ResidueView& residue_view) {
    const Molecule& source = *residue_view.molecule;
    const Residue& source_residue = residue_view.get();
    const AtomId atom_begin = source_residue.atom_begin;
    const AtomId atom_end = atom_begin + source_residue.atom_count;
    const auto contains = [atom_begin, atom_end](AtomId atom) {
        return atom >= atom_begin && atom < atom_end;
    };
    const auto all_contained = [&contains](std::initializer_list<AtomId> atoms) {
        return std::all_of(atoms.begin(), atoms.end(), contains);
    };
    const auto local = [atom_begin](AtomId atom) { return atom - atom_begin; };

    auto molecule = std::make_shared<Molecule>(source_residue.name);
    Residue residue = source_residue;
    residue.atom_begin = 0;
    molecule->residues.push_back(std::move(residue));
    molecule->atoms.assign(source.atoms.begin() + atom_begin, source.atoms.begin() + atom_end);
    for (auto& atom : molecule->atoms) {
        atom.residue = 0;
    }

    for (const auto& bond : source.explicit_bonds) {
        if (all_contained({bond.atom1, bond.atom2})) {
            molecule->explicit_bonds.push_back({local(bond.atom1), local(bond.atom2)});
        }
    }
    for (const auto& link : source.residue_links) {
        if (all_contained({link.atom1, link.atom2})) {
            molecule->residue_links.push_back({local(link.atom1), local(link.atom2)});
        }
    }
    for (const auto& item : source.virtual_atoms) {
        if (all_contained({item.virtual_atom, item.atom0, item.atom1, item.atom2})) {
            molecule->virtual_atoms.push_back(
                {local(item.virtual_atom), local(item.atom0), local(item.atom1), local(item.atom2),
                 item.k1, item.k2});
        }
    }
    for (const auto& item : source.harmonic_impropers) {
        if (all_contained({item.atom0, item.atom1, item.atom2, item.atom3})) {
            molecule->harmonic_impropers.push_back(
                {local(item.atom0), local(item.atom1), local(item.atom2), local(item.atom3),
                 item.k, item.phi0});
        }
    }
    molecule->cmap_types = source.cmap_types;
    for (const auto& item : source.cmaps) {
        if (all_contained({item.atom0, item.atom1, item.atom2, item.atom3, item.atom4})) {
            molecule->cmaps.push_back(
                {local(item.atom0), local(item.atom1), local(item.atom2), local(item.atom3),
                 local(item.atom4), item.type});
        }
    }
    for (const auto& item : source.nb14_extras) {
        if (all_contained({item.atom1, item.atom2})) {
            molecule->nb14_extras.push_back(
                {local(item.atom1), local(item.atom2), item.a, item.b, item.kee});
        }
    }
    for (const auto& item : source.urey_bradleys) {
        if (all_contained({item.atom0, item.atom1, item.atom2})) {
            molecule->urey_bradleys.push_back(
                {local(item.atom0), local(item.atom1), local(item.atom2),
                 item.k, item.b, item.k_ub, item.r13});
        }
    }
    for (const auto& item : source.ryckaert_bellemans) {
        if (all_contained({item.atom0, item.atom1, item.atom2, item.atom3})) {
            molecule->ryckaert_bellemans.push_back(
                {local(item.atom0), local(item.atom1), local(item.atom2), local(item.atom3),
                 item.c0, item.c1, item.c2, item.c3, item.c4, item.c5});
        }
    }
    for (const auto& item : source.soft_bonds) {
        if (all_contained({item.atom1, item.atom2})) {
            molecule->soft_bonds.push_back(
                {local(item.atom1), local(item.atom2), item.k, item.b, item.from_a_or_b});
        }
    }

    if (source.topology_override) {
        Topology topology;
        for (const auto& item : source.topology_override->bonds) {
            if (all_contained({item.atom1, item.atom2})) {
                topology.bonds.push_back(
                    {local(item.atom1), local(item.atom2), item.k, item.length});
            }
        }
        for (const auto& item : source.topology_override->angles) {
            if (all_contained({item.atom1, item.atom2, item.atom3})) {
                topology.angles.push_back(
                    {local(item.atom1), local(item.atom2), local(item.atom3), item.k, item.theta});
            }
        }
        for (const auto& item : source.topology_override->dihedrals) {
            if (all_contained({item.atom1, item.atom2, item.atom3, item.atom4})) {
                topology.dihedrals.push_back(
                    {local(item.atom1), local(item.atom2), local(item.atom3), local(item.atom4),
                     item.periodicity, item.k, item.phase});
            }
        }
        topology.exclusions.resize(source_residue.atom_count);
        for (AtomId atom = atom_begin; atom < atom_end; ++atom) {
            if (atom >= source.topology_override->exclusions.size()) {
                break;
            }
            for (const AtomId excluded : source.topology_override->exclusions[atom]) {
                if (contains(excluded)) {
                    topology.exclusions[local(atom)].push_back(local(excluded));
                }
            }
        }
        for (const auto& item : source.topology_override->nb14s) {
            if (all_contained({item.atom1, item.atom2})) {
                topology.nb14s.push_back(
                    {local(item.atom1), local(item.atom2), item.k_lj, item.k_ee});
            }
        }
        molecule->topology_override = std::move(topology);
    }

    molecule->listed_force_definitions = source.listed_force_definitions;
    molecule->sw_parameters = source.sw_parameters;
    molecule->edip_parameters = source.edip_parameters;
    molecule->box_length = source.box_length;
    molecule->box_angle = source.box_angle;
    molecule->has_box = source.has_box;
    molecule->has_gb_parameters = source.has_gb_parameters;
    molecule->write_min_bonded_parameters = source.write_min_bonded_parameters;
    molecule->write_subsys_division = source.write_subsys_division;
    molecule->write_lj_soft_core = source.write_lj_soft_core;
    molecule->ignore_missing_atoms = source.ignore_missing_atoms;
    if (!molecule->validate()) {
        throw std::runtime_error("internal error: invalid molecule materialized from residue");
    }
    return molecule;
}

std::shared_ptr<Molecule> save_residue_bundle_object(
    const ResidueView& residue, py::object prefix, py::object dirname, py::object protocol) {
    require_empty_native_protocol(protocol);
    auto molecule = molecule_from_residue_view(residue);
    require_native_bundle_compatibility(molecule);
    save_sponge_input_bundle(*molecule, normalize_bundle_prefix(prefix, molecule->name),
                             py::str(dirname).cast<std::string>());
    return molecule;
}

std::shared_ptr<Molecule> save_template_like_bundle_object(
    py::object source, py::object prefix, py::object dirname, py::object protocol) {
    require_empty_native_protocol(protocol);
    if (!py::hasattr(source, "name")) {
        throw py::type_error("save_sponge_input_bundle expects a Molecule or residue template");
    }
    auto molecule = std::make_shared<Molecule>(
        get_template_molecule(py::str(source.attr("name")).cast<std::string>()));
    require_native_bundle_compatibility(molecule);
    save_sponge_input_bundle(*molecule, normalize_bundle_prefix(prefix, molecule->name),
                             py::str(dirname).cast<std::string>());
    return molecule;
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

void bind_core_module(py::module_& m) {
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
        .def("add_bond_soft", &Molecule::add_bond_soft, py::arg("atom1"), py::arg("atom2"), py::arg("k"),
             py::arg("b"), py::arg("from_AorB"))
        .def("Add_Bond_Soft", &Molecule::add_bond_soft, py::arg("atom1"), py::arg("atom2"), py::arg("k"),
             py::arg("b"), py::arg("from_AorB"))
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
        .def("set_ignore_missing_atoms", &Molecule::set_ignore_missing_atoms, py::arg("enabled") = true)
        .def("Set_Ignore_Missing_Atoms", &Molecule::set_ignore_missing_atoms, py::arg("enabled") = true)
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

    m.def("add_molecule", &add_molecule_object, py::arg("molecule"), py::arg("other"));
    m.def("replace_residues", &replace_residues_object, py::arg("molecule"), py::arg("replacements"),
          py::arg("residue_sort_keys") = std::vector<double>{}, py::arg("sort") = true);
    m.def("reorder_atoms_by_template", &reorder_atoms_by_template_object, py::arg("molecule"),
          py::arg("template_molecule"));
    m.def("set_box_padding", &set_box_padding_object, py::arg("molecule"), py::arg("padding") = 0.5,
          py::arg("center") = true);
    m.def("save_sponge_input", &save_sponge_input_object, py::arg("molecule"), py::arg("prefix") = "",
          py::arg("dirname") = ".");
    m.def("save_sponge_input_bundle", &save_sponge_input_bundle_object, py::arg("molecule"),
          py::arg("prefix") = py::none(), py::arg("dirname") = ".", py::kw_only(),
          py::arg("protocol") = py::none());
    m.def("save_sponge_input_bundle", &save_residue_bundle_object, py::arg("molecule"),
          py::arg("prefix") = py::none(), py::arg("dirname") = ".", py::kw_only(),
          py::arg("protocol") = py::none());
    m.def("save_sponge_input_bundle", &save_residuetype_bundle_object, py::arg("molecule"),
          py::arg("prefix") = py::none(), py::arg("dirname") = ".", py::kw_only(),
          py::arg("protocol") = py::none());
    m.def("save_sponge_input_bundle", &save_template_like_bundle_object, py::arg("molecule"),
          py::arg("prefix") = py::none(), py::arg("dirname") = ".", py::kw_only(),
          py::arg("protocol") = py::none());
    m.def("merge_dual_topology", &merge_dual_topology_object, py::arg("molecule"), py::arg("residue_index"),
          py::arg("residue_b_molecule"), py::arg("match_b_to_a"));
    m.def("merge_force_field", &merge_force_field_object, py::arg("molecule_a"), py::arg("molecule_b"),
          py::arg("default_lambda"), py::arg("specific_lambda") = std::unordered_map<std::string, double>{});
    m.def("get_template_molecule", &get_template_molecule_object);
    m.def("molecule_from_residuetype", [](const ResidueType& residue_type) {
        auto molecule = std::make_shared<Molecule>(residue_type.name());
        molecule->append_residue_from_type(residue_type, 0.0, 0.0, 0.0);
        return molecule;
    }, py::arg("residue_type"));
}

}  // namespace xpongecpp
