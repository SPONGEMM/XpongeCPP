#pragma once

#include "core.hpp"

#include <fstream>
#include <memory>
#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace xpongecpp {

inline std::string read_python_input(py::object source) {
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

inline std::vector<ResidueView> residue_views(const std::shared_ptr<Molecule>& molecule) {
    std::vector<ResidueView> views;
    views.reserve(molecule->residue_count());
    for (ResidueId i = 0; i < molecule->residue_count(); ++i) {
        views.push_back({molecule, i});
    }
    return views;
}

inline std::vector<AtomView> residue_atom_views(const ResidueView& residue_view) {
    std::vector<AtomView> views;
    const auto& residue = residue_view.get();
    views.reserve(residue.atom_count);
    for (std::uint32_t local = 0; local < residue.atom_count; ++local) {
        views.push_back({residue_view.molecule, residue.atom_begin + local});
    }
    return views;
}

inline std::vector<AtomView> molecule_atom_views(const std::shared_ptr<Molecule>& molecule) {
    std::vector<AtomView> views;
    views.reserve(molecule->atom_count());
    for (AtomId i = 0; i < molecule->atom_count(); ++i) {
        views.push_back({molecule, i});
    }
    return views;
}

inline std::vector<std::array<AtomId, 2>> molecule_explicit_bonds(const std::shared_ptr<Molecule>& molecule) {
    std::vector<std::array<AtomId, 2>> bonds;
    bonds.reserve(molecule->explicit_bonds.size());
    for (const auto& bond : molecule->explicit_bonds) {
        bonds.push_back({bond.atom1, bond.atom2});
    }
    return bonds;
}

inline std::vector<std::array<AtomId, 2>> molecule_residue_links(const std::shared_ptr<Molecule>& molecule) {
    std::vector<std::array<AtomId, 2>> links;
    links.reserve(molecule->residue_links.size());
    for (const auto& link : molecule->residue_links) {
        links.push_back({link.atom1, link.atom2});
    }
    return links;
}

void bind_core_module(py::module_& m);
void bind_io_module(py::module_& m);
void bind_forcefield_module(py::module_& m);
void bind_assign_module(py::module_& m);

}  // namespace xpongecpp
