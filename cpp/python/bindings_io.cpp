#include "bindings_internal.hpp"
#include "pdb_internal.hpp"

namespace xpongecpp {
namespace {

char selector_char(py::object value, char fallback) {
    const auto text = pdb_trimmed_copy(py::cast<std::string>(py::str(value)));
    return text.empty() ? fallback : text[0];
}

py::object dict_get_first(const py::dict& dict, std::initializer_list<const char*> keys) {
    for (const char* key : keys) {
        if (dict.contains(py::str(key))) {
            return dict[py::str(key)];
        }
    }
    return py::none();
}

bool dict_get_bool(const py::dict& dict, std::initializer_list<const char*> keys) {
    py::object value = dict_get_first(dict, keys);
    return !value.is_none() && py::cast<bool>(value);
}

template <typename Options>
void add_terminal_selector(Options& options, py::object chain, py::object resseq, py::object insertion,
                           bool n_terminal, bool c_terminal) {
    if (resseq.is_none()) {
        throw std::invalid_argument("terminal residue selector requires a residue sequence number");
    }
    PdbLoadOptions::TerminalResidue selector;
    selector.chain_id = selector_char(chain, ' ');
    selector.resseq = py::cast<int>(resseq);
    selector.insertion_code = selector_char(insertion, ' ');
    selector.n_terminal = n_terminal;
    selector.c_terminal = c_terminal;
    options.terminal_residues.push_back(selector);
}

template <typename Options>
void add_terminal_kind_selector(Options& options, py::object chain, py::object resseq, py::object insertion,
                                py::object terminal_kind) {
    auto kind = pdb_upper_copy(pdb_trimmed_copy(py::cast<std::string>(py::str(terminal_kind))));
    if (kind == "N" || kind == "HEAD" || kind == "N_TERMINAL" || kind == "N-TERMINAL") {
        add_terminal_selector(options, chain, resseq, insertion, true, false);
    } else if (kind == "C" || kind == "TAIL" || kind == "C_TERMINAL" || kind == "C-TERMINAL") {
        add_terminal_selector(options, chain, resseq, insertion, false, true);
    } else {
        throw std::invalid_argument("invalid terminal kind: " + kind);
    }
}

template <typename Options>
void parse_terminal_residues(Options& options, py::object terminal_residues) {
    if (terminal_residues.is_none()) {
        return;
    }
    for (py::handle item_handle : terminal_residues) {
        py::object item = py::reinterpret_borrow<py::object>(item_handle);
        if (py::isinstance<py::dict>(item)) {
            py::dict dict = py::cast<py::dict>(item);
            py::object chain = dict_get_first(dict, {"chain_id", "chainId"});
            py::object resseq = dict_get_first(dict, {"res_seq", "resSeq", "residue_seq", "residueSeq"});
            py::object insertion = dict_get_first(dict, {"icode", "insertion_code", "insertionCode"});
            if (chain.is_none()) {
                chain = py::str("");
            }
            if (insertion.is_none()) {
                insertion = py::str("");
            }
            const bool n_terminal = dict_get_bool(dict, {"n_terminal", "nTerminal"});
            const bool c_terminal = dict_get_bool(dict, {"c_terminal", "cTerminal"});
            if (n_terminal || c_terminal) {
                add_terminal_selector(options, chain, resseq, insertion, n_terminal, c_terminal);
            }
            py::object kind = dict_get_first(dict, {"terminal", "kind", "place"});
            if (!kind.is_none() && !pdb_trimmed_copy(py::cast<std::string>(py::str(kind))).empty()) {
                add_terminal_kind_selector(options, chain, resseq, insertion, kind);
            }
            continue;
        }
        if (py::isinstance<py::str>(item) || !py::isinstance<py::sequence>(item)) {
            throw std::invalid_argument("invalid terminal residue selector");
        }
        py::sequence sequence = py::cast<py::sequence>(item);
        if (sequence.size() == 4) {
            add_terminal_kind_selector(options, sequence[0], sequence[1], sequence[2], sequence[3]);
            continue;
        }
        if (sequence.size() == 5) {
            add_terminal_selector(options, sequence[0], sequence[1], sequence[2], py::cast<bool>(sequence[3]),
                                  py::cast<bool>(sequence[4]));
            continue;
        }
        throw std::invalid_argument("invalid terminal residue selector");
    }
}

MmcifResidueLinkAtom parse_mmcif_residue_link_atom(py::object item) {
    if (!py::isinstance<py::dict>(item)) {
        throw std::invalid_argument("mmCIF residue link atom selector must be a dict");
    }
    py::dict dict = py::cast<py::dict>(item);
    py::object chain = dict_get_first(dict, {"chain_id", "chainId"});
    py::object resseq = dict_get_first(dict, {"residue_seq", "residueSeq", "res_seq", "resSeq"});
    py::object insertion = dict_get_first(dict, {"insertion_code", "insertionCode", "icode"});
    py::object residue_name = dict_get_first(dict, {"residue_name", "residueName", "resname"});
    py::object atom_name = dict_get_first(dict, {"atom_name", "atomName", "name"});
    if (chain.is_none()) chain = py::str("");
    if (insertion.is_none()) insertion = py::str("");
    if (resseq.is_none() || atom_name.is_none()) {
        throw std::invalid_argument("mmCIF residue link atom selector requires residue_seq and atom_name");
    }
    MmcifResidueLinkAtom atom;
    atom.chain_id = selector_char(chain, ' ');
    atom.resseq = py::cast<int>(resseq);
    atom.insertion_code = selector_char(insertion, ' ');
    atom.residue_name = residue_name.is_none() ? std::string{} : py::cast<std::string>(py::str(residue_name));
    atom.atom_name = py::cast<std::string>(py::str(atom_name));
    return atom;
}

void parse_mmcif_residue_links(MmcifLoadOptions& options, py::object residue_links) {
    if (residue_links.is_none()) {
        return;
    }
    for (py::handle item_handle : residue_links) {
        py::object item = py::reinterpret_borrow<py::object>(item_handle);
        if (!py::isinstance<py::dict>(item)) {
            throw std::invalid_argument("mmCIF residue link selector must be a dict");
        }
        py::dict dict = py::cast<py::dict>(item);
        py::object atom1 = dict_get_first(dict, {"atom_a", "atomA", "atom1"});
        py::object atom2 = dict_get_first(dict, {"atom_b", "atomB", "atom2"});
        if (atom1.is_none() || atom2.is_none()) {
            throw std::invalid_argument("mmCIF residue link selector requires atom_a and atom_b");
        }
        options.residue_links.push_back({parse_mmcif_residue_link_atom(atom1), parse_mmcif_residue_link_atom(atom2)});
    }
}

std::shared_ptr<Molecule> load_pdb_object(py::object source, bool judge_histone, const std::string& position_need,
                                          bool ignore_hydrogen, bool ignore_unknown_name, bool ignore_seqres,
                                          bool ignore_conect, bool read_cryst1, py::object unterminal_residues,
                                          py::object terminal_residues, bool infer_terminals) {
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
    parse_terminal_residues(options, terminal_residues);
    options.infer_terminals = infer_terminals;
    return std::make_shared<Molecule>(load_pdb_text(read_python_input(source), options));
}

std::shared_ptr<Molecule> load_mmcif_object(py::object source, bool judge_histone, const std::string& position_need,
                                            bool ignore_hydrogen, bool ignore_unknown_name, bool ignore_seqres,
                                            bool read_cell, py::object unterminal_residues,
                                            py::object terminal_residues, bool infer_terminals, py::object model_id,
                                            py::object residue_links) {
    MmcifLoadOptions options;
    options.judge_histone = judge_histone;
    options.position_need = position_need.empty() ? 'A' : position_need[0];
    options.ignore_hydrogen = ignore_hydrogen;
    options.ignore_unknown_name = ignore_unknown_name;
    options.ignore_seqres = ignore_seqres;
    options.read_cell = read_cell;
    if (!unterminal_residues.is_none()) {
        for (const auto item : unterminal_residues) {
            options.unterminal_residues.push_back(py::str(item));
        }
    }
    parse_terminal_residues(options, terminal_residues);
    options.infer_terminals = infer_terminals;
    if (!model_id.is_none()) {
        options.model_id = py::cast<std::string>(py::str(model_id));
    }
    parse_mmcif_residue_links(options, residue_links);
    return std::make_shared<Molecule>(load_mmcif_text(read_python_input(source), options));
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
          py::arg("unterminal_residues") = py::none(), py::arg("terminal_residues") = py::none(),
          py::arg("infer_terminals") = true);
    m.def("load_mmcif", &load_mmcif_object, py::arg("source"), py::arg("judge_histone") = true,
          py::arg("position_need") = "A", py::arg("ignore_hydrogen") = false,
          py::arg("ignore_unknown_name") = false, py::arg("ignore_seqres") = true,
          py::arg("read_cell") = true, py::arg("unterminal_residues") = py::none(),
          py::arg("terminal_residues") = py::none(), py::arg("infer_terminals") = true,
          py::arg("model_id") = py::none(), py::arg("residue_links") = py::none());
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
