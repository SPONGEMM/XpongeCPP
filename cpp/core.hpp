#pragma once

#include <array>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace xpongecpp {

using AtomId = std::uint32_t;
using ResidueId = std::uint32_t;

struct Atom {
    std::string name;
    std::string type;
    std::string element;
    ResidueId residue{0};
    std::int32_t serial{0};
    char altloc{' '};
    double occupancy{1.0};
    double temp_factor{0.0};
    std::string record_name{"ATOM"};
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double charge{0.0};
    double mass{0.0};
    std::string lj_type_b;
    std::string sw_type;
    std::string edip_type;
    double gb_radius{0.0};
    double gb_scaler{0.0};
    int subsys{0};
    bool bad_coordinate{false};
    bool zero_lj_atom{false};
};

struct Residue {
    std::string name;
    std::string type_name;
    std::string original_name;
    char chain_id{' '};
    char effective_chain_id{' '};
    std::uint32_t segment_id{0};
    std::int32_t pdb_resseq{0};
    char insertion_code{' '};
    bool is_hetero{false};
    AtomId atom_begin{0};
    std::uint32_t atom_count{0};
};

struct ResidueLink {
    AtomId atom1{0};
    AtomId atom2{0};
};

struct Bond {
    AtomId atom1{0};
    AtomId atom2{0};
    double k{0.0};
    double length{0.0};
};

struct Angle {
    AtomId atom1{0};
    AtomId atom2{0};
    AtomId atom3{0};
    double k{0.0};
    double theta{0.0};
};

struct Dihedral {
    AtomId atom1{0};
    AtomId atom2{0};
    AtomId atom3{0};
    AtomId atom4{0};
    int periodicity{1};
    double k{0.0};
    double phase{0.0};
};

struct NB14 {
    AtomId atom1{0};
    AtomId atom2{0};
    double k_lj{0.5};
    double k_ee{0.833333};
};

struct VirtualAtom2 {
    AtomId virtual_atom{0};
    AtomId atom0{0};
    AtomId atom1{0};
    AtomId atom2{0};
    double k1{0.0};
    double k2{0.0};
};

struct HarmonicImproper {
    AtomId atom0{0};
    AtomId atom1{0};
    AtomId atom2{0};
    AtomId atom3{0};
    double k{0.0};
    double phi0{0.0};
};

struct CMapType {
    std::uint32_t resolution{0};
    std::vector<double> parameters;
};

struct CMap {
    AtomId atom0{0};
    AtomId atom1{0};
    AtomId atom2{0};
    AtomId atom3{0};
    AtomId atom4{0};
    std::uint32_t type{0};
};

struct NB14Extra {
    AtomId atom1{0};
    AtomId atom2{0};
    double a{0.0};
    double b{0.0};
    double kee{0.0};
};

struct UreyBradley {
    AtomId atom0{0};
    AtomId atom1{0};
    AtomId atom2{0};
    double k{0.0};
    double b{0.0};
    double k_ub{0.0};
    double r13{0.0};
};

struct RyckaertBellemans {
    AtomId atom0{0};
    AtomId atom1{0};
    AtomId atom2{0};
    AtomId atom3{0};
    double c0{0.0};
    double c1{0.0};
    double c2{0.0};
    double c3{0.0};
    double c4{0.0};
    double c5{0.0};
};

struct SoftBond {
    AtomId atom1{0};
    AtomId atom2{0};
    double k{0.0};
    double b{0.0};
    int from_a_or_b{0};
};

struct StillingerWeberParameter {
    double a_big{0.0};
    double b_big{0.0};
    double epsilon{0.0};
    double p{0.0};
    double q{0.0};
    double a{0.0};
    double gamma{0.0};
    double sigma{0.0};
    double lambda{0.0};
    double b{0.0};
};

struct EDIPParameter {
    double a_big{0.0};
    double b_big{0.0};
    double a{0.0};
    double c{0.0};
    double alpha{0.0};
    double beta{0.0};
    double eta{0.0};
    double gamma{0.0};
    double lambda{0.0};
    double mu{0.0};
    double rho{0.0};
    double sigma{0.0};
    double q0{0.0};
    double u1{0.0};
    double u2{0.0};
    double u3{0.0};
    double u4{0.0};
};

struct DihedralTerm {
    int periodicity{1};
    double k{0.0};
    double phase{0.0};
};

struct NB14Scale {
    double k_lj{0.5};
    double k_ee{0.833333};
};

struct BondTerm {
    double k{0.0};
    double length{0.0};
};

struct AngleTerm {
    double k{0.0};
    double theta{0.0};
};

struct AmberImproperMatch {
    DihedralTerm term;
    bool exact{false};
    int wildcard_count{0};
};

enum class LJCombiningRule {
    LorentzBerthelot,
    GoodHope,
};

struct Topology {
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;
    std::vector<std::vector<AtomId>> exclusions;
    std::vector<NB14> nb14s;
};

struct ResidueTypeAtom {
    std::string name;
    std::string type;
    std::string element;
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double charge{0.0};
    double mass{0.0};
};

struct ResidueTypeBond {
    std::uint32_t atom1{0};
    std::uint32_t atom2{0};
};

class ResidueType {
public:
    explicit ResidueType(std::string name);

    const std::string& name() const noexcept;
    std::uint64_t version() const noexcept;
    std::size_t atom_count() const noexcept;
    std::size_t bond_count() const noexcept;
    const std::vector<ResidueTypeAtom>& atoms() const noexcept;
    const std::vector<ResidueTypeBond>& bonds() const noexcept;
    const std::string& head() const noexcept;
    const std::string& tail() const noexcept;
    const std::string& head_next() const noexcept;
    const std::string& tail_next() const noexcept;
    double head_length() const noexcept;
    double tail_length() const noexcept;
    const std::unordered_map<std::string, std::string>& connect_atoms() const noexcept;

    void add_atom(const std::string& name, const std::string& type, double x, double y, double z,
                  double charge, double mass);
    void add_connectivity(const std::string& atom1, const std::string& atom2);
    std::uint32_t atom_index(const std::string& name) const;
    void set_head(const std::string& atom, double length = 1.5, const std::string& next = "");
    void set_tail(const std::string& atom, double length = 1.5, const std::string& next = "");
    void set_connect_atom(const std::string& key, const std::string& atom);

private:
    std::string name_;
    std::uint64_t version_{0};
    std::vector<ResidueTypeAtom> atoms_;
    std::vector<ResidueTypeBond> bonds_;
    std::unordered_map<std::string, std::uint32_t> atom_name_to_index_;
    std::string head_;
    std::string tail_;
    std::string head_next_;
    std::string tail_next_;
    double head_length_{1.5};
    double tail_length_{1.5};
    std::unordered_map<std::string, std::string> connect_atoms_;
};

struct PdbLoadOptions {
    bool judge_histone{true};
    char position_need{'A'};
    bool ignore_hydrogen{false};
    bool ignore_unknown_name{false};
    bool ignore_seqres{true};
    bool ignore_conect{true};
    bool read_cryst1{true};
    std::vector<std::string> unterminal_residues;
    struct TerminalResidue {
        char chain_id{' '};
        int resseq{0};
        char insertion_code{' '};
        bool n_terminal{false};
        bool c_terminal{false};
    };
    std::vector<TerminalResidue> terminal_residues;
    bool infer_terminals{true};
};

struct MmcifResidueLinkAtom {
    char chain_id{' '};
    int resseq{0};
    char insertion_code{' '};
    std::string residue_name;
    std::string atom_name;
};

struct MmcifResidueLink {
    MmcifResidueLinkAtom atom1;
    MmcifResidueLinkAtom atom2;
};

struct MmcifLoadOptions {
    bool judge_histone{true};
    char position_need{'A'};
    bool ignore_hydrogen{false};
    bool ignore_unknown_name{false};
    bool ignore_seqres{true};
    bool read_cell{true};
    std::vector<std::string> unterminal_residues;
    std::vector<PdbLoadOptions::TerminalResidue> terminal_residues;
    bool infer_terminals{true};
    std::optional<std::string> model_id;
    std::vector<MmcifResidueLink> residue_links;
};

class Molecule {
public:
    explicit Molecule(std::string name = "MOL");

    std::string name;
    std::vector<Atom> atoms;
    std::vector<Residue> residues;
    std::vector<ResidueLink> explicit_bonds;
    std::vector<ResidueLink> residue_links;
    std::vector<VirtualAtom2> virtual_atoms;
    std::vector<HarmonicImproper> harmonic_impropers;
    std::vector<CMapType> cmap_types;
    std::vector<CMap> cmaps;
    std::vector<NB14Extra> nb14_extras;
    std::vector<UreyBradley> urey_bradleys;
    std::vector<RyckaertBellemans> ryckaert_bellemans;
    std::vector<SoftBond> soft_bonds;
    std::vector<std::string> listed_force_definitions;
    std::unordered_map<std::string, StillingerWeberParameter> sw_parameters;
    std::unordered_map<std::string, EDIPParameter> edip_parameters;
    std::optional<Topology> topology_override;
    std::array<double, 3> box_length{0.0, 0.0, 0.0};
    std::array<double, 3> box_angle{90.0, 90.0, 90.0};
    bool has_box{false};
    bool has_gb_parameters{false};
    bool write_min_bonded_parameters{false};
    bool write_subsys_division{false};
    bool write_lj_soft_core{false};
    bool ignore_missing_atoms{true};

    std::size_t atom_count() const noexcept;
    std::size_t residue_count() const noexcept;
    const Atom& atom(AtomId id) const;
    Atom& atom(AtomId id);
    const Residue& residue(ResidueId id) const;
    Residue& residue(ResidueId id);
    void append_residue_from_type(const ResidueType& type, double dx, double dy, double dz);
    void add_molecule(const Molecule& other);
    void add_molecule_linked(const Molecule& other, bool link);
    void add_residue_link(AtomId atom1, AtomId atom2);
    void add_virtual_atom2(AtomId virtual_atom, AtomId atom0, AtomId atom1, AtomId atom2, double k1, double k2);
    void add_improper_dihedral(AtomId atom0, AtomId atom1, AtomId atom2, AtomId atom3, double k, double phi0);
    std::uint32_t add_cmap_type(std::uint32_t resolution, const std::vector<double>& parameters);
    void add_cmap(AtomId atom0, AtomId atom1, AtomId atom2, AtomId atom3, AtomId atom4, std::uint32_t type);
    void add_nb14_extra(AtomId atom1, AtomId atom2, double a, double b, double kee);
    void add_urey_bradley(AtomId atom0, AtomId atom1, AtomId atom2, double k, double b, double k_ub, double r13);
    void add_ryckaert_bellemans(AtomId atom0, AtomId atom1, AtomId atom2, AtomId atom3,
                                double c0, double c1, double c2, double c3, double c4, double c5);
    void add_bond_soft(AtomId atom1, AtomId atom2, double k, double b, int from_a_or_b);
    void add_listed_force_definition(const std::string& definition);
    void set_gb_radius(const std::string& radius_set = "modified_bondi_radii");
    void enable_min_bonded_parameters(bool enabled = true) noexcept;
    void enable_subsys_division(bool enabled = true) noexcept;
    void enable_lj_soft_core(bool enabled = true) noexcept;
    void add_sw_type(const std::string& name, double a_big, double b_big, double epsilon, double p, double q,
                     double a, double gamma, double sigma, double lambda, double b);
    void add_edip_type(const std::string& name, double a_big, double b_big, double a, double c, double alpha,
                       double beta, double eta, double gamma, double lambda, double mu, double rho, double sigma,
                       double q0, double u1, double u2, double u3, double u4);
    void set_ignore_missing_atoms(bool enabled = true) noexcept;
    void set_box_padding(double padding, bool center);
    void replace_residues(const std::unordered_map<ResidueId, Molecule>& replacements,
                          const std::vector<double>& residue_sort_keys = {}, bool sort = true);
    void reorder_atoms_by_template(const Molecule& template_molecule);
    bool validate() const;
    std::unordered_map<std::string, std::size_t> residue_counts() const;
};

class Assign {
public:
    explicit Assign(std::string name = "ASN");

    std::string name;
    std::vector<std::string> elements;
    std::vector<std::string> element_details;
    std::vector<std::string> names;
    std::vector<std::array<double, 3>> coordinates;
    std::vector<double> charges;
    std::vector<int> formal_charges;
    std::vector<std::unordered_map<std::uint32_t, int>> bonds;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> bond_sequence;
    std::vector<std::unordered_map<std::uint32_t, std::set<std::string>>> bond_markers;
    std::vector<std::unordered_map<std::string, int>> atom_markers;
    std::vector<std::string> atom_types;
    std::vector<std::vector<std::uint32_t>> rings;
    bool built{false};

    void add_atom(const std::string& element, double x, double y, double z,
                  const std::string& name = "", double charge = 0.0);
    void add_bond(std::uint32_t atom1, std::uint32_t atom2, int order);
    void delete_bond(std::uint32_t atom1, std::uint32_t atom2);
    void delete_atom(std::uint32_t atom);
    void set_charge(std::uint32_t atom, double charge);
    void set_charges(const std::vector<double>& new_charges);
    void set_formal_charge(std::uint32_t atom, int charge);
    void set_coordinate(std::uint32_t atom, double x, double y, double z);
    void set_atom_type(std::uint32_t atom, const std::string& atom_type);
    void determine_connectivity(double simple_cutoff);
    bool check_connectivity() const;
    bool atom_judge(std::uint32_t atom, const std::string& mask) const;
    bool atom_judge(std::uint32_t atom, const std::vector<std::string>& masks) const;
    void add_atom_marker(std::uint32_t atom, const std::string& marker);
    void add_bond_marker(std::uint32_t atom1, std::uint32_t atom2, const std::string& marker, bool only1 = false);
    bool has_atom_marker(std::uint32_t atom, const std::string& marker) const;
    int atom_marker_count(std::uint32_t atom, const std::string& marker) const;
    bool has_bond_marker(std::uint32_t atom1, std::uint32_t atom2, const std::string& marker) const;
    bool determine_bond_order(bool check_formal_charge = true, std::optional<int> total_charge = std::nullopt);
    bool determine_bond_order_custom(
        bool check_formal_charge,
        std::optional<int> total_charge,
        int max_step,
        int max_stat,
        const std::vector<std::vector<std::pair<int, int>>>& penalty_scores,
        const std::function<bool(const Assign&)>& extra_criteria = {});
    void determine_ring_and_bond_type();
    void kekulize();
    void determine_atom_type(const std::string& rule);
    void calculate_tpacm4_charge(const std::string& atom_type_table,
                                 const std::string& charge_table,
                                 int total_charge);
    ResidueType to_residuetype(const std::string& name) const;
    Molecule to_molecule(const std::string& name) const;
    std::size_t atom_count() const noexcept;
    std::size_t bond_count() const noexcept;
};

Molecule load_pdb_text(const std::string& text);
Molecule load_pdb_text(const std::string& text, const PdbLoadOptions& options);
Molecule load_mmcif_text(const std::string& text);
Molecule load_mmcif_text(const std::string& text, const MmcifLoadOptions& options);
Molecule load_mol2_text(const std::string& text);
Assign get_assignment_from_mol2_text(const std::string& text,
                                     std::optional<int> total_charge = std::nullopt,
                                     bool total_charge_from_partial_sum = false);
Assign get_assignment_from_xyz_text(const std::string& text);
Assign get_assignment_from_pdb_text(const std::string& text);
Assign get_assignment_from_residuetype(const ResidueType& residue_type);
std::string assignment_to_mol2_text(const Assign& assignment, const std::string& residue_name);
std::string assignment_to_pdb_text(const Assign& assignment, const std::string& residue_name);
std::vector<std::string> implemented_gaff_assign_types();
std::vector<std::string> implemented_gaff2_assign_types();
std::vector<std::array<double, 3>> generate_resp_mk_grid(
    const std::vector<std::string>& atoms,
    const std::vector<std::array<double, 3>>& atom_coordinates_bohr,
    double area_density = 1.0,
    int layer = 4,
    const std::unordered_map<std::string, double>& radius = {}
);
struct RespFitDebugResult {
    std::vector<double> esp_charges;
    std::vector<double> stage1_charges;
    std::vector<double> final_charges;
    std::vector<int> stage2_restrained_groups;
    std::unordered_map<std::string, double> timings;
};
std::vector<double> fit_resp_from_esp_cpp(
    const Assign& assign,
    const std::vector<std::array<double, 3>>& atom_coordinates_bohr,
    const std::vector<double>& nuclear_charges,
    const std::vector<std::array<double, 3>>& grid_points_bohr,
    const std::vector<double>& esp_values_au,
    int charge,
    const std::vector<std::vector<int>>& extra_equivalence = {},
    double a1 = 0.0005,
    double a2 = 0.001,
    bool two_stage = true,
    bool only_esp = false
);
RespFitDebugResult fit_resp_from_esp_cpp_debug(
    const Assign& assign,
    const std::vector<std::array<double, 3>>& atom_coordinates_bohr,
    const std::vector<double>& nuclear_charges,
    const std::vector<std::array<double, 3>>& grid_points_bohr,
    const std::vector<double>& esp_values_au,
    int charge,
    const std::vector<std::vector<int>>& extra_equivalence = {},
    double a1 = 0.0005,
    double a2 = 0.001,
    bool two_stage = true,
    bool only_esp = false
);
void add_solvent_box(Molecule& molecule, const Molecule& solvent, double distance, double tolerance,
                     std::int64_t n_solvent, std::uint64_t seed = 0);
void add_solvent_box(Molecule& molecule, const Molecule& solvent, const std::array<double, 6>& distance,
                     double tolerance, std::int64_t n_solvent, std::uint64_t seed = 0);
void add_ions(Molecule& molecule, const std::unordered_map<std::string, std::int64_t>& counts,
              std::uint64_t seed = 0, const std::string& solvent_residue = "WAT");
std::unordered_map<std::string, std::filesystem::path> save_sponge_input(Molecule& molecule,
                                                                         const std::string& prefix,
                                                                         const std::filesystem::path& dirname);
void save_pdb(const Molecule& molecule, const std::filesystem::path& filename);
void save_mol2(const Molecule& molecule, const std::filesystem::path& filename);
struct GroData {
    std::vector<std::array<double, 3>> coordinates;
    std::array<double, 6> box{0.0, 0.0, 0.0, 90.0, 90.0, 90.0};
};
struct CoordinateData {
    std::vector<std::array<double, 3>> coordinates;
    std::array<double, 6> box{0.0, 0.0, 0.0, 90.0, 90.0, 90.0};
    bool has_box{false};
};
struct PsfData {
    Molecule molecule;
    std::unordered_map<std::string, Molecule> molecules;
};
CoordinateData load_coordinate_text(const std::string& text);
CoordinateData load_coordinate_text(const std::string& text, Molecule& molecule);
CoordinateData load_rst7_text(const std::string& text);
CoordinateData load_rst7_text(const std::string& text, Molecule& molecule);
GroData load_gro_text(const std::string& text);
GroData load_gro_text(const std::string& text, bool read_box_angle);
GroData load_gro_text(const std::string& text, Molecule& molecule, bool read_box_angle = true);
void save_gro(const Molecule& molecule, const std::filesystem::path& filename);
PsfData load_molpsf_text(const std::string& text, const std::string& split_by = "connectivity");
Molecule load_gromacs_topology_file(const std::filesystem::path& filename);
Molecule load_opls_itp_file(const std::filesystem::path& filename);
void load_charmm_parameter_file(const std::filesystem::path& filename);
Molecule load_charmm_topology_file(const std::filesystem::path& filename);
void load_sw_parameter_file(const std::filesystem::path& filename, Molecule& molecule);
void load_edip_parameter_file(const std::filesystem::path& filename, Molecule& molecule);

Topology build_topology(const Molecule& molecule);

void register_ff14sb();
void register_tip3p();
void register_amber_parmdat_file(const std::filesystem::path& filename);
void register_amber_frcmod_file(const std::filesystem::path& filename);
void register_amber_lj_parameter(const std::string& atom_type, const std::string& lj_type, double epsilon, double rmin);
void register_amber_bond_parameter(const std::string& atom_type1, const std::string& atom_type2, double k,
                                   double length);
void register_amber_angle_parameter(const std::array<std::string, 3>& atom_types, double k, double theta);
void register_amber_proper_dihedral_parameter(const std::array<std::string, 4>& atom_types, int periodicity,
                                              double k, double phase, bool reset = false);
void register_amber_improper_dihedral_parameter(const std::array<std::string, 4>& atom_types, int periodicity,
                                                double k, double phase);
void register_amber_nb14_scale(const std::string& atom_type1, const std::string& atom_type4, double k_lj,
                               double k_ee);
void register_amber_cmap_parameter(const std::string& key, std::uint32_t resolution,
                                   const std::vector<double>& parameters);
void clear_amber_dihedral_parameters();
void clear_amber_improper_parameters();
void set_lj_combining_rule(LJCombiningRule rule);
LJCombiningRule lj_combining_rule();
bool has_amber_cmap_parameters();
void apply_amber_cmaps(Molecule& molecule);
std::vector<DihedralTerm> find_amber_proper_terms(const std::array<std::string, 4>& atom_types);
std::optional<DihedralTerm> find_amber_improper_term(const std::array<std::string, 4>& atom_types);
std::optional<AmberImproperMatch> find_amber_improper_match(const std::array<std::string, 4>& atom_types);
std::optional<NB14Scale> find_amber_nb14_scale(const std::string& atom_type1, const std::string& atom_type4);
std::optional<NB14Scale> find_amber_nb14_dihedral_scale(
    const std::array<std::string, 4>& atom_types);
std::optional<BondTerm> find_amber_bond_term(const std::string& atom_type1, const std::string& atom_type2);
std::optional<AngleTerm> find_amber_angle_term(const std::array<std::string, 3>& atom_types);
std::string find_amber_lj_type(const std::string& atom_type);
std::optional<std::pair<double, double>> find_amber_lj_parameter(const std::string& lj_type);
std::optional<double> find_amber_atom_type_mass(const std::string& atom_type);
std::optional<double> find_external_atom_type_mass(const std::string& atom_type);
std::pair<Molecule, Molecule> merge_dual_topology(const Molecule& molecule, ResidueId residue_index,
                                                  const Molecule& residue_b_molecule,
                                                  const std::unordered_map<std::uint32_t, std::uint32_t>& match_b_to_a);
Molecule merge_force_field(const Molecule& molecule_a, const Molecule& molecule_b, double default_lambda,
                           const std::unordered_map<std::string, double>& specific_lambda);
void register_residue_templates_from_mol2_text(const std::string& text);
void register_residue_templates_from_mol2_file(const std::filesystem::path& filename);
void register_template_molecule_from_mol2_file(const std::filesystem::path& filename);
void register_template_virtual_atom2(const std::string& template_name, const std::string& virtual_atom,
                                     const std::string& atom0, const std::string& atom1, const std::string& atom2,
                                     double k1, double k2);
void configure_residue_template_head(const std::string& template_name, const std::string& atom,
                                     double length = 1.5, const std::string& next = "");
void configure_residue_template_tail(const std::string& template_name, const std::string& atom,
                                     double length = 1.5, const std::string& next = "");
void configure_residue_template_connect_atom(const std::string& template_name, const std::string& key,
                                             const std::string& atom);
void register_residue_template_alias(const std::string& alias_name, const std::string& template_name);
bool has_template(const std::string& name);
std::size_t template_atom_count(const std::string& name);
std::vector<std::string> registered_template_names();
Molecule get_template_molecule(const std::string& name);
const ResidueType& get_residue_template(const std::string& name);
void register_pdb_residue_name_mapping(const std::string& place, const std::string& pdb_name,
                                       const std::string& real_name);
void register_pdb_residue_alias_mapping(const std::string& pdb_name, const std::string& real_name);
void register_his_mapping(const std::string& residue_name, const std::string& hid, const std::string& hie,
                          const std::string& hip);

double default_mass_for_element(const std::string& element);
std::string guess_element_from_mass(double mass);
std::string guess_element(const std::string& atom_name, const std::string& explicit_element = "");

}  // namespace xpongecpp
