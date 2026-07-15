#include "core.hpp"
#include "sponge_writers.hpp"

#include <hdf5.h>
#include <highfive/highfive.hpp>

#include <chrono>
#include <cmath>
#include <stdexcept>

namespace xpongecpp {
namespace {

class H5File {
public:
  explicit H5File(const std::filesystem::path &path)
      : handle(path.string(), HighFive::File::Overwrite) {}
  H5File(const H5File &) = delete;
  H5File &operator=(const H5File &) = delete;
  HighFive::File handle;
};

void ensure_groups(H5File &file, const std::string &dataset_path) {
  std::size_t offset = 1;
  while (true) {
    const auto slash = dataset_path.find('/', offset);
    if (slash == std::string::npos)
      break;
    const auto group = dataset_path.substr(0, slash);
    if (!file.handle.exist(group))
      file.handle.createGroup(group);
    offset = slash + 1;
  }
}

template <typename T>
void write_array(H5File &file, const std::string &path,
                 const std::vector<hsize_t> &dimensions,
                 const std::vector<T> &values) {
  ensure_groups(file, path);
  auto dataset =
      file.handle.createDataSet<T>(path, HighFive::DataSpace(dimensions));
  if (!values.empty())
    dataset.write_raw(values.data());
}

template <typename T>
void write_scalar(H5File &file, const std::string &path, T value) {
  ensure_groups(file, path);
  auto dataset =
      file.handle.createDataSet<T>(path, HighFive::DataSpace::From(value));
  dataset.write(value);
}

void write_string(H5File &file, const std::string &path,
                  const std::string &value) {
  ensure_groups(file, path);
  auto dataset = file.handle.createDataSet<std::string>(
      path, HighFive::DataSpace::From(value));
  dataset.write(value);
}

void write_strings(H5File &file, const std::string &path,
                   const std::vector<std::string> &values) {
  ensure_groups(file, path);
  auto dataset = file.handle.createDataSet<std::string>(
      path, HighFive::DataSpace({values.size()}));
  if (!values.empty())
    dataset.write(values);
}

void create_hard_link(H5File &file, const std::string &target,
                      const std::string &link_path) {
  ensure_groups(file, link_path);
  if (file.handle.exist(link_path))
    return;
  if (H5Lcreate_hard(file.handle.getId(), target.c_str(), file.handle.getId(),
                     link_path.c_str(), H5P_DEFAULT, H5P_DEFAULT) < 0) {
    throw std::runtime_error("failed to create HDF5 hard link: " + link_path);
  }
}

template <typename T>
void set_attribute(H5File &file, const std::string &object_path,
                   const std::string &name, const T &value) {
  auto attribute =
      file.handle.getDataSet(object_path)
          .createAttribute<T>(name, HighFive::DataSpace::From(value));
  attribute.write(value);
}

template <typename T>
void set_group_attribute(H5File &file, const std::string &object_path,
                         const std::string &name, const T &value) {
  auto attribute =
      file.handle.getGroup(object_path)
          .createAttribute<T>(name, HighFive::DataSpace::From(value));
  attribute.write(value);
}

template <typename T>
void set_group_array_attribute(H5File &file, const std::string &object_path,
                               const std::string &name,
                               const std::vector<T> &values) {
  auto attribute =
      file.handle.getGroup(object_path)
          .createAttribute<T>(name, HighFive::DataSpace({values.size()}));
  attribute.write(values);
}

void finalize_topology(H5File &file, std::size_t atom_count) {
  const std::string version = "xponge.legacy_to_bundle.v1";
  const std::string hash = "xpongecpp-native-v1:" + std::to_string(atom_count);
  write_string(file, "/schema/name", "sponge.topology.h5");
  write_string(file, "/schema/version", version);
  write_string(file, "/parameters/sponge/schema/name", "sponge.topology.h5");
  write_string(file, "/parameters/sponge/schema/version", version);
  write_string(file, "/topology/atom_order_hash", hash);
  write_string(file, "/topology/topology_hash", hash);
  write_string(file, "/topology/forcefield_hash", hash);
}

std::pair<float, float> lj_ab(const std::string &lhs, const std::string &rhs) {
  const auto a = find_amber_lj_parameter(lhs);
  const auto b = find_amber_lj_parameter(rhs);
  if (!a || !b)
    throw std::runtime_error("missing Amber LJ parameter for bundled topology");
  const double epsilon = std::sqrt(a->first * b->first);
  double radius = a->second + b->second;
  if (lj_combining_rule() == LJCombiningRule::GoodHope)
    radius = 2.0 * std::sqrt(a->second * b->second);
  const double r6 = std::pow(radius, 6.0);
  return {static_cast<float>(epsilon * r6 * r6),
          static_cast<float>(2.0 * epsilon * r6)};
}

void write_native_topology(H5File &file, const Molecule &molecule) {
  const auto topology = build_topology(molecule);
  const auto atom_count = molecule.atoms.size();
  std::vector<float> mass, charge;
  std::vector<std::int32_t> residue_index;
  std::vector<std::string> atom_name, atom_type_name, residue_name;
  for (const auto &atom : molecule.atoms) {
    mass.push_back(static_cast<float>(atom.mass));
    charge.push_back(static_cast<float>(atom.charge * 18.2223));
    residue_index.push_back(static_cast<std::int32_t>(atom.residue));
    atom_name.push_back(atom.name);
    atom_type_name.push_back(atom.type);
  }
  std::vector<std::int64_t> residue_offset{0};
  for (const auto &residue : molecule.residues) {
    residue_offset.push_back(residue_offset.back() + residue.atom_count);
    residue_name.push_back(residue.name);
  }
  write_array(file, "/atoms/mass", {atom_count}, mass);
  write_array(file, "/atoms/charge", {atom_count}, charge);
  set_attribute(file, "/atoms/charge", "unit", std::string("Amber"));
  write_array(file, "/atoms/residue_index", {atom_count}, residue_index);
  write_array(file, "/residues/atom_offset", {residue_offset.size()},
              residue_offset);
  write_strings(file, "/parameters/xponge/atoms/name", atom_name);
  write_strings(file, "/parameters/xponge/atoms/type_name", atom_type_name);
  write_strings(file, "/parameters/xponge/residues/name", residue_name);

  std::vector<std::int32_t> bond_atoms, angle_atoms, dihedral_atoms,
      dihedral_periodicity, nb14_atoms;
  std::vector<float> bond_k, bond_r0, angle_k, angle_theta0, dihedral_k,
      dihedral_phi0, nb14_params;
  for (const auto &item : topology.bonds) {
    bond_atoms.insert(bond_atoms.end(),
                      {static_cast<std::int32_t>(item.atom1),
                       static_cast<std::int32_t>(item.atom2)});
    bond_k.push_back(item.k);
    bond_r0.push_back(item.length);
  }
  for (const auto &item : topology.angles) {
    angle_atoms.insert(angle_atoms.end(),
                       {static_cast<std::int32_t>(item.atom1),
                        static_cast<std::int32_t>(item.atom2),
                        static_cast<std::int32_t>(item.atom3)});
    angle_k.push_back(item.k);
    angle_theta0.push_back(item.theta);
  }
  for (const auto &item : topology.dihedrals) {
    dihedral_atoms.insert(dihedral_atoms.end(),
                          {static_cast<std::int32_t>(item.atom1),
                           static_cast<std::int32_t>(item.atom2),
                           static_cast<std::int32_t>(item.atom3),
                           static_cast<std::int32_t>(item.atom4)});
    dihedral_periodicity.push_back(item.periodicity);
    dihedral_k.push_back(item.k);
    dihedral_phi0.push_back(item.phase);
  }
  for (const auto &item : topology.nb14s) {
    nb14_atoms.insert(nb14_atoms.end(),
                      {static_cast<std::int32_t>(item.atom1),
                       static_cast<std::int32_t>(item.atom2)});
    if (molecule.write_lj_soft_core) {
      const auto lhs = find_amber_lj_type(molecule.atoms[item.atom1].type);
      const auto rhs = find_amber_lj_type(molecule.atoms[item.atom2].type);
      const auto [a, b] = lj_ab(lhs, rhs);
      nb14_params.insert(nb14_params.end(),
                         {static_cast<float>(item.k_lj * a * 12.0),
                          static_cast<float>(item.k_lj * b * 6.0),
                          static_cast<float>(item.k_ee)});
    } else {
      nb14_params.insert(nb14_params.end(), {static_cast<float>(item.k_lj),
                                             static_cast<float>(item.k_ee)});
    }
  }
  write_array(file, "/forcefield/bond/atoms", {topology.bonds.size(), 2},
              bond_atoms);
  write_array(file, "/forcefield/bond/k", {topology.bonds.size()}, bond_k);
  write_array(file, "/forcefield/bond/r0", {topology.bonds.size()}, bond_r0);
  write_scalar<std::int64_t>(file, "/forcefield/bond/count",
                             topology.bonds.size());
  write_array(file, "/forcefield/angle/atoms", {topology.angles.size(), 3},
              angle_atoms);
  write_array(file, "/forcefield/angle/k", {topology.angles.size()}, angle_k);
  write_array(file, "/forcefield/angle/theta0", {topology.angles.size()},
              angle_theta0);
  write_scalar<std::int64_t>(file, "/forcefield/angle/count",
                             topology.angles.size());
  write_array(file, "/forcefield/dihedral/atoms",
              {topology.dihedrals.size(), 4}, dihedral_atoms);
  write_array(file, "/forcefield/dihedral/periodicity",
              {topology.dihedrals.size()}, dihedral_periodicity);
  write_array(file, "/forcefield/dihedral/k", {topology.dihedrals.size()},
              dihedral_k);
  write_array(file, "/forcefield/dihedral/phi0", {topology.dihedrals.size()},
              dihedral_phi0);
  write_scalar<std::int64_t>(file, "/forcefield/dihedral/count",
                             topology.dihedrals.size());
  std::vector<std::int32_t> improper_atoms;
  std::vector<float> improper_pk, improper_phi0;
  for (const auto &item : molecule.harmonic_impropers) {
    improper_atoms.insert(improper_atoms.end(),
                          {static_cast<std::int32_t>(item.atom2),
                           static_cast<std::int32_t>(item.atom0),
                           static_cast<std::int32_t>(item.atom1),
                           static_cast<std::int32_t>(item.atom3)});
    improper_pk.push_back(item.k);
    improper_phi0.push_back(item.phi0);
  }
  if (!molecule.harmonic_impropers.empty()) {
    write_array(file, "/forcefield/improper/atoms",
                {molecule.harmonic_impropers.size(), 4}, improper_atoms);
    write_array(file, "/forcefield/improper/pk",
                {molecule.harmonic_impropers.size()}, improper_pk);
    write_array(file, "/forcefield/improper/phi0",
                {molecule.harmonic_impropers.size()}, improper_phi0);
    write_scalar<std::int64_t>(file, "/forcefield/improper/count",
                               molecule.harmonic_impropers.size());
  }
  write_array(file, "/forcefield/nb14/atoms", {topology.nb14s.size(), 2},
              nb14_atoms);
  write_array(file, "/forcefield/nb14/params",
              {topology.nb14s.size(), molecule.write_lj_soft_core ? 3UL : 2UL},
              nb14_params);
  write_scalar<std::int64_t>(file, "/forcefield/nb14/count",
                             topology.nb14s.size());

  std::vector<std::int32_t> nb14_extra_atoms;
  std::vector<float> nb14_extra_params;
  for (const auto &item : molecule.nb14_extras) {
    nb14_extra_atoms.insert(nb14_extra_atoms.end(),
                            {static_cast<std::int32_t>(item.atom1),
                             static_cast<std::int32_t>(item.atom2)});
    nb14_extra_params.insert(nb14_extra_params.end(),
                             {static_cast<float>(item.a),
                              static_cast<float>(item.b),
                              static_cast<float>(item.kee)});
  }
  if (!molecule.nb14_extras.empty()) {
    write_array(file, "/forcefield/nb14_extra/atoms",
                {molecule.nb14_extras.size(), 2}, nb14_extra_atoms);
    write_array(file, "/forcefield/nb14_extra/params",
                {molecule.nb14_extras.size(), 3}, nb14_extra_params);
  }

  if (!molecule.virtual_atoms.empty()) {
    std::vector<std::int32_t> types, atoms, from;
    std::vector<std::int64_t> from_offset{0}, parameter_offset{0};
    std::vector<float> parameters;
    for (const auto &item : molecule.virtual_atoms) {
      types.push_back(2);
      atoms.push_back(static_cast<std::int32_t>(item.virtual_atom));
      from.insert(from.end(), {static_cast<std::int32_t>(item.atom0),
                               static_cast<std::int32_t>(item.atom1),
                               static_cast<std::int32_t>(item.atom2)});
      parameters.insert(parameters.end(), {static_cast<float>(item.k1),
                                           static_cast<float>(item.k2)});
      from_offset.push_back(from.size());
      parameter_offset.push_back(parameters.size());
    }
    write_array(file, "/forcefield/virtual_atom/type", {types.size()}, types);
    write_array(file, "/forcefield/virtual_atom/atom", {atoms.size()}, atoms);
    write_array(file, "/forcefield/virtual_atom/from_offset",
                {from_offset.size()}, from_offset);
    write_array(file, "/forcefield/virtual_atom/from", {from.size()}, from);
    write_array(file, "/forcefield/virtual_atom/parameter_offset",
                {parameter_offset.size()}, parameter_offset);
    write_array(file, "/forcefield/virtual_atom/parameter", {parameters.size()},
                parameters);
    write_scalar<std::int64_t>(file, "/forcefield/virtual_atom/count",
                               types.size());
  }

  if (!molecule.urey_bradleys.empty()) {
    std::vector<std::int32_t> atoms;
    std::vector<float> angle_k, angle_theta0, bond_k, bond_r0;
    for (const auto &item : molecule.urey_bradleys) {
      atoms.insert(atoms.end(), {static_cast<std::int32_t>(item.atom0),
                                 static_cast<std::int32_t>(item.atom1),
                                 static_cast<std::int32_t>(item.atom2)});
      angle_k.push_back(item.k);
      angle_theta0.push_back(item.b);
      bond_k.push_back(item.k_ub);
      bond_r0.push_back(item.r13);
    }
    write_array(file, "/forcefield/urey_bradley/atoms",
                {molecule.urey_bradleys.size(), 3}, atoms);
    write_array(file, "/forcefield/urey_bradley/angle_k", {angle_k.size()},
                angle_k);
    write_array(file, "/forcefield/urey_bradley/angle_theta0",
                {angle_theta0.size()}, angle_theta0);
    write_array(file, "/forcefield/urey_bradley/bond_k", {bond_k.size()},
                bond_k);
    write_array(file, "/forcefield/urey_bradley/bond_r0", {bond_r0.size()},
                bond_r0);
  }

  if (molecule.has_gb_parameters) {
    std::vector<float> gb_params;
    gb_params.reserve(atom_count * 2);
    for (const auto &atom : molecule.atoms) {
      gb_params.insert(gb_params.end(), {static_cast<float>(atom.gb_radius),
                                         static_cast<float>(atom.gb_scaler)});
    }
    write_array(file, "/forcefield/gb/params", {atom_count, 2}, gb_params);
  }
  if (molecule.write_subsys_division) {
    std::vector<std::int32_t> subsystem;
    subsystem.reserve(atom_count);
    for (const auto &atom : molecule.atoms)
      subsystem.push_back(atom.subsys);
    write_array(file, "/forcefield/subsys_division", {atom_count}, subsystem);
  }

  if (!molecule.cmaps.empty()) {
    std::vector<std::int32_t> atoms, type, resolution;
    std::vector<float> grid_value;
    std::vector<std::uint32_t> source_types;
    std::unordered_map<std::uint32_t, std::int32_t> output_type;
    for (const auto &item : molecule.cmaps) {
      auto [it, inserted] = output_type.emplace(item.type, output_type.size());
      if (inserted)
        source_types.push_back(item.type);
      atoms.insert(atoms.end(), {static_cast<std::int32_t>(item.atom0),
                                 static_cast<std::int32_t>(item.atom1),
                                 static_cast<std::int32_t>(item.atom2),
                                 static_cast<std::int32_t>(item.atom3),
                                 static_cast<std::int32_t>(item.atom4)});
      type.push_back(it->second);
    }
    for (const auto source_type : source_types) {
      if (source_type >= molecule.cmap_types.size()) {
        throw std::runtime_error("CMAP references an undefined type");
      }
      const auto &item = molecule.cmap_types[source_type];
      resolution.push_back(item.resolution);
      for (const auto value : item.parameters)
        grid_value.push_back(value);
    }
    write_array(file, "/forcefield/cmap/atoms", {molecule.cmaps.size(), 5},
                atoms);
    write_array(file, "/forcefield/cmap/type", {type.size()}, type);
    write_array(file, "/forcefield/cmap/resolution", {resolution.size()},
                resolution);
    write_array(file, "/forcefield/cmap/grid_value", {grid_value.size()},
                grid_value);
    write_scalar<std::int64_t>(file, "/forcefield/cmap/count",
                               molecule.cmaps.size());
  }

  std::vector<std::int64_t> exclude_offset{0};
  std::vector<std::int32_t> exclude_list;
  for (std::size_t atom_index = 0; atom_index < topology.exclusions.size();
       ++atom_index) {
    for (const auto excluded_atom : topology.exclusions[atom_index]) {
      // SPONGE's legacy exclusion payload stores each unordered pair once,
      // under the lower-index atom.  build_topology() intentionally exposes a
      // symmetric adjacency list, so materialize only its upper triangle.
      if (atom_index < excluded_atom)
        exclude_list.push_back(excluded_atom);
    }
    exclude_offset.push_back(exclude_list.size());
  }
  write_array(file, "/topology/exclusions/offset", {exclude_offset.size()},
              exclude_offset);
  write_array(file, "/topology/exclusions/list", {exclude_list.size()},
              exclude_list);

  auto collect_lj_state = [&](bool state_b) {
    std::vector<std::string> names;
    std::unordered_map<std::string, std::int32_t> index;
    std::vector<std::int32_t> atom_type;
    for (const auto &atom : molecule.atoms) {
      const auto source_type =
          state_b && !atom.lj_type_b.empty() ? atom.lj_type_b : atom.type;
      const auto name = find_amber_lj_type(source_type);
      if (!index.count(name)) {
        index[name] = static_cast<std::int32_t>(names.size());
        names.push_back(name);
      }
      atom_type.push_back(index.at(name));
    }
    std::vector<float> pair_a, pair_b;
    for (std::size_t i = 0; i < names.size(); ++i)
      for (std::size_t j = 0; j <= i; ++j) {
        const auto [a, b] = lj_ab(names[i], names[j]);
        pair_a.push_back(a);
        pair_b.push_back(b);
      }
    return std::make_tuple(std::move(names), std::move(atom_type),
                           std::move(pair_a), std::move(pair_b));
  };
  auto [lj_names, atom_lj_type, pair_a, pair_b] = collect_lj_state(false);
  if (molecule.write_lj_soft_core) {
    auto [lj_names_b, atom_lj_type_b, pair_b_a, pair_b_b] =
        collect_lj_state(true);
    write_scalar<std::int32_t>(file,
                               "/forcefield/lj_soft_core/atom_type_count_A",
                               lj_names.size());
    write_scalar<std::int32_t>(file,
                               "/forcefield/lj_soft_core/atom_type_count_B",
                               lj_names_b.size());
    write_array(file, "/forcefield/lj_soft_core/pair_AA", {pair_a.size()},
                pair_a);
    write_array(file, "/forcefield/lj_soft_core/pair_AB", {pair_b.size()},
                pair_b);
    write_array(file, "/forcefield/lj_soft_core/pair_BA", {pair_b_a.size()},
                pair_b_a);
    write_array(file, "/forcefield/lj_soft_core/pair_BB", {pair_b_b.size()},
                pair_b_b);
    write_array(file, "/forcefield/lj_soft_core/atom_type_A", {atom_count},
                atom_lj_type);
    write_array(file, "/forcefield/lj_soft_core/atom_type_B", {atom_count},
                atom_lj_type_b);
  } else {
    write_scalar<std::int32_t>(file, "/forcefield/lj/atom_type_count",
                               lj_names.size());
    write_array(file, "/forcefield/lj/pair_A_12", {pair_a.size()}, pair_a);
    write_array(file, "/forcefield/lj/pair_B_6", {pair_b.size()}, pair_b);
    write_array(file, "/forcefield/lj/type", {atom_count}, atom_lj_type);
  }
  write_scalar<std::int64_t>(file, "/topology/atom_count", atom_count);
}

void write_protocol(const std::filesystem::path &path, std::size_t atom_count) {
  H5File file(path);
  const std::string version = "xponge.legacy_to_bundle.v1";
  const std::string hash = "xpongecpp-native-v1:" + std::to_string(atom_count);
  write_string(file, "/schema/name", "sponge.protocol.h5");
  write_string(file, "/schema/version", version);
  write_string(file, "/parameters/sponge/schema/name", "sponge.protocol.h5");
  write_string(file, "/parameters/sponge/schema/version", version);
  write_string(file, "/protocol/topology_compatibility/topology_hash", hash);
  write_string(file, "/identity/content_hash", "xpongecpp-native-protocol-v1");
  write_scalar<std::int64_t>(file, "/protocol/cv_count", 0);
  write_scalar<std::int64_t>(file, "/protocol/restraint_count", 0);
}

void write_restart(const Molecule &molecule,
                   const std::filesystem::path &path) {
  H5File file(path);
  for (const auto &group :
       {"/h5md", "/h5md/creator", "/run", "/particles/all",
        "/parameters/restart", "/parameters/restart/rng_state",
        "/parameters/restart/integrator_state",
        "/parameters/restart/thermostat", "/parameters/restart/barostat",
        "/parameters/restart/protocol_sidecars", "/parameters/restart/bias",
        "/parameters/restart/bias/sits", "/parameters/restart/bias/meta"}) {
    if (!file.handle.exist(group))
      file.handle.createGroup(group);
  }
  const auto shift = sponge_coordinate_shift_for_export(molecule);
  std::vector<float> positions;
  positions.reserve(molecule.atoms.size() * 3);
  for (const auto &atom : molecule.atoms) {
    positions.push_back(static_cast<float>(atom.x + shift[0]));
    positions.push_back(static_cast<float>(atom.y + shift[1]));
    positions.push_back(static_cast<float>(atom.z + shift[2]));
  }
  const auto box = sponge_coordinate_box_for_export(molecule);
  const double alpha = box[3] * 3.14159265358979323846 / 180.0;
  const double beta = box[4] * 3.14159265358979323846 / 180.0;
  const double gamma = box[5] * 3.14159265358979323846 / 180.0;
  std::vector<float> edges{
      static_cast<float>(box[0]),
      0,
      0,
      static_cast<float>(box[1] * std::cos(gamma)),
      static_cast<float>(box[1] * std::sin(gamma)),
      0,
      static_cast<float>(box[2] * std::cos(beta)),
      static_cast<float>(box[2] *
                         (std::cos(alpha) - std::cos(beta) * std::cos(gamma)) /
                         std::sin(gamma)),
      0};
  const double z2 = box[2] * box[2] - edges[6] * edges[6] - edges[7] * edges[7];
  edges[8] = static_cast<float>(std::sqrt(std::max(0.0, z2)));
  write_array(file, "/particles/all/position/value",
              {1, molecule.atoms.size(), 3}, positions);
  write_array(file, "/particles/all/box/edges/value", {1, 3, 3}, edges);
  write_array<std::int64_t>(file, "/particles/all/step", {1}, {0});
  write_array<double>(file, "/particles/all/time", {1}, {0.0});
  set_group_array_attribute<std::int32_t>(file, "/h5md", "version", {1, 1});
  set_group_attribute(file, "/h5md/creator", "name", std::string("XpongeCPP"));
  set_group_attribute(file, "/h5md/creator", "version", std::string("0.1.2"));
  set_attribute(file, "/particles/all/time", "unit", std::string("ps"));
  set_attribute(file, "/particles/all/position/value", "unit",
                std::string("Angstrom"));
  set_attribute(file, "/particles/all/box/edges/value", "unit",
                std::string("Angstrom"));
  set_group_attribute<std::int32_t>(file, "/particles/all/box", "dimension", 3);
  set_group_array_attribute<std::string>(file, "/particles/all/box", "boundary",
                                         {"periodic", "periodic", "periodic"});
  create_hard_link(file, "/particles/all/step", "/particles/all/position/step");
  create_hard_link(file, "/particles/all/time", "/particles/all/position/time");
  create_hard_link(file, "/particles/all/step",
                   "/particles/all/box/edges/step");
  create_hard_link(file, "/particles/all/time",
                   "/particles/all/box/edges/time");
  write_string(file, "/parameters/sponge/schema/name", "sponge.restart.h5");
  write_string(file, "/parameters/sponge/schema/version",
               "xponge.legacy_to_bundle.v1");
  write_string(file, "/parameters/sponge/output/status", "finalized");
  write_scalar<std::int64_t>(file, "/parameters/sponge/output/frame_count", 1);
  write_scalar<std::int64_t>(file,
                             "/parameters/sponge/output/last_complete_step", 0);
  write_array<double>(file, "/parameters/sponge/output/last_complete_time", {1},
                      {0.0});
  write_scalar<std::int64_t>(file, "/run/current_step", 0);
  write_array<double>(file, "/run/current_time", {1}, {0.0});
  write_strings(file, "/parameters/sponge/output/particle_streams", {"all"});
  write_string(file, "/run/state_type", "restart");
}

class TemporaryBundleFiles {
public:
  TemporaryBundleFiles(const std::filesystem::path &topology_target,
                       const std::filesystem::path &protocol_target,
                       const std::filesystem::path &restart_target)
      : targets{topology_target, protocol_target, restart_target} {
    const auto nonce =
        std::chrono::steady_clock::now().time_since_epoch().count();
    for (const auto &target : targets) {
      temporary.emplace_back(target.string() + ".tmp." + std::to_string(nonce));
    }
  }

  ~TemporaryBundleFiles() {
    if (committed)
      return;
    for (const auto &path : temporary) {
      std::error_code ignored;
      std::filesystem::remove(path, ignored);
    }
  }

  void commit() {
    for (std::size_t i = 0; i < targets.size(); ++i) {
      std::error_code error;
      std::filesystem::rename(temporary[i], targets[i], error);
      if (error) {
        std::filesystem::remove(targets[i], error);
        error.clear();
        std::filesystem::rename(temporary[i], targets[i], error);
      }
      if (error) {
        throw std::runtime_error("failed to install bundled input file: " +
                                 targets[i].string() + ": " + error.message());
      }
    }
    committed = true;
  }

  std::array<std::filesystem::path, 3> targets;
  std::vector<std::filesystem::path> temporary;
  bool committed{false};
};

} // namespace

std::unordered_map<std::string, std::filesystem::path>
save_sponge_input_bundle(const Molecule &input_molecule,
                         const std::string &prefix,
                         const std::filesystem::path &dirname) {
  std::optional<Molecule> molecule_with_generated_cmaps;
  if (input_molecule.cmaps.empty() && has_amber_cmap_parameters()) {
    molecule_with_generated_cmaps = input_molecule;
    apply_amber_cmaps(*molecule_with_generated_cmaps);
  }
  const Molecule &molecule = molecule_with_generated_cmaps
                                 ? *molecule_with_generated_cmaps
                                 : input_molecule;
  if (!molecule.validate())
    throw std::invalid_argument("cannot export invalid molecule bundle");
  const std::filesystem::path relative(prefix.empty() ? molecule.name : prefix);
  if (relative.is_absolute())
    throw std::invalid_argument("bundle prefix must be relative");
  std::filesystem::create_directories(dirname);
  const auto base = std::filesystem::weakly_canonical(dirname / relative);
  const auto root = std::filesystem::weakly_canonical(dirname);
  const auto relative_to_root = base.lexically_relative(root);
  if (relative_to_root.empty() || *relative_to_root.begin() == "..") {
    throw std::invalid_argument("bundle prefix escapes output directory");
  }
  const auto topology_path =
      std::filesystem::path(base.string() + "_topology.spgt.h5");
  const auto protocol_path =
      std::filesystem::path(base.string() + "_protocol.spgp.h5");
  const auto restart_path =
      std::filesystem::path(base.string() + "_restart.spgr.h5");
  std::filesystem::create_directories(topology_path.parent_path());

  TemporaryBundleFiles files(topology_path, protocol_path, restart_path);

  {
    H5File topology(files.temporary[0]);
    write_native_topology(topology, molecule);
    finalize_topology(topology, molecule.atoms.size());
  }
  write_protocol(files.temporary[1], molecule.atoms.size());
  write_restart(molecule, files.temporary[2]);
  files.commit();
  return {{"topology", topology_path},
          {"protocol", protocol_path},
          {"restart", restart_path}};
}

} // namespace xpongecpp
