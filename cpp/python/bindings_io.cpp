#include "bindings_internal.hpp"

namespace xpongecpp {
namespace {

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

}  // namespace

void bind_io_module(py::module_& m) {
    m.def("load_pdb", &load_pdb_object, py::arg("source"), py::arg("judge_histone") = true,
          py::arg("position_need") = "A", py::arg("ignore_hydrogen") = false,
          py::arg("ignore_unknown_name") = false, py::arg("ignore_seqres") = true,
          py::arg("ignore_conect") = true, py::arg("read_cryst1") = true,
          py::arg("unterminal_residues") = py::none());
    m.def("load_mol2", &load_mol2_object);
    m.def("load_molpsf", &load_molpsf_object, py::arg("source"), py::arg("split_by") = "connectivity");
    m.def("save_pdb", &save_pdb_object, py::arg("molecule"), py::arg("filename"), py::arg("write_cryst1") = true);
    m.def("save_mol2", &save_mol2_object, py::arg("molecule"), py::arg("filename"));
    m.def("load_coordinate", &load_coordinate_object, py::arg("source"), py::arg("molecule") = py::none());
    m.def("load_rst7", &load_rst7_object, py::arg("source"), py::arg("molecule") = py::none());
    m.def("load_gro", &load_gro_object, py::arg("source"), py::arg("molecule") = py::none(),
          py::arg("read_box_angle") = true);
    m.def("save_gro", &save_gro_object, py::arg("molecule"), py::arg("filename"));
}

}  // namespace xpongecpp
