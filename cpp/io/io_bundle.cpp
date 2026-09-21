#include "core.hpp"
#include "sha256.hpp"
#include "sponge_writers.hpp"

#include <hdf5.h>
#include <highfive/highfive.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <type_traits>

namespace xpongecpp {
namespace {

constexpr const char *kInputSchemaVersion = "sponge.input.v2";
constexpr const char *kRyckaertBellemansParameterText =
    "int atom_a, int atom_b, int atom_c, int atom_d, float c0, float c1, "
    "float c2, float c3, float c4, float c5";
constexpr const char *kRyckaertBellemansPotential =
    "SADfloat<15> cphi = cosf(phi_abcd - CONSTANT_Pi);\n"
    "SADfloat<15> cphi2 = cphi * cphi;\n"
    "SADfloat<15> cphi3 = cphi2 * cphi;\n"
    "SADfloat<15> cphi4 = cphi3 * cphi;\n"
    "SADfloat<15> cphi5 = cphi4 * cphi;\n"
    "E = c0 + c1 * cphi + c2 * cphi2 + c3 * cphi3 + c4 * cphi4 + "
    "c5 * cphi5;";

std::string generate_uuid_v4() {
  std::random_device source;
  std::array<std::uint8_t, 16> bytes{};
  for (auto &byte : bytes)
    byte = static_cast<std::uint8_t>(source());
  bytes[6] = static_cast<std::uint8_t>((bytes[6] & 0x0fU) | 0x40U);
  bytes[8] = static_cast<std::uint8_t>((bytes[8] & 0x3fU) | 0x80U);
  std::ostringstream value;
  value << std::hex << std::setfill('0');
  for (std::size_t index = 0; index < bytes.size(); ++index) {
    if (index == 4 || index == 6 || index == 8 || index == 10)
      value << '-';
    value << std::setw(2) << static_cast<unsigned int>(bytes[index]);
  }
  return value.str();
}

struct TrackedDataset {
  std::string dtype;
  std::vector<std::size_t> dimensions;
  std::vector<std::uint8_t> bytes;
};

template <typename T> const char *numpy_dtype_name();
template <> const char *numpy_dtype_name<float>() { return "float32"; }
template <> const char *numpy_dtype_name<double>() { return "float64"; }
template <> const char *numpy_dtype_name<std::int32_t>() { return "int32"; }
template <> const char *numpy_dtype_name<std::int64_t>() { return "int64"; }
template <> const char *numpy_dtype_name<std::uint32_t>() { return "uint32"; }

std::string numpy_shape(const std::vector<std::size_t> &dimensions) {
  if (dimensions.empty())
    return "()";
  std::string output = "(";
  for (std::size_t index = 0; index < dimensions.size(); ++index) {
    if (index)
      output += ", ";
    output += std::to_string(dimensions[index]);
  }
  if (dimensions.size() == 1)
    output += ',';
  return output + ')';
}

class DatasetHashTracker {
public:
  template <typename T>
  void add_numeric(const std::string &path,
                   const std::vector<std::size_t> &dimensions,
                   const std::vector<T> &values) {
    static_assert(std::is_arithmetic_v<T>, "numeric dataset required");
    TrackedDataset dataset;
    dataset.dtype = numpy_dtype_name<T>();
    dataset.dimensions = dimensions;
    dataset.bytes.resize(values.size() * sizeof(T));
    if (!values.empty())
      std::memcpy(dataset.bytes.data(), values.data(), dataset.bytes.size());
    datasets_[path] = std::move(dataset);
  }

  template <typename T>
  void add_scalar(const std::string &path, T value) {
    add_numeric(path, {}, std::vector<T>{value});
  }

  void add_strings(const std::string &path,
                   const std::vector<std::string> &values) {
    TrackedDataset dataset;
    dataset.dtype = "object";
    dataset.dimensions = {values.size()};
    std::size_t byte_count = values.empty() ? 0 : values.size() - 1;
    for (const auto &value : values)
      byte_count += value.size();
    dataset.bytes.reserve(byte_count);
    for (std::size_t index = 0; index < values.size(); ++index) {
      if (index)
        dataset.bytes.push_back(0);
      dataset.bytes.insert(dataset.bytes.end(), values[index].begin(),
                           values[index].end());
    }
    datasets_[path] = std::move(dataset);
  }

  void add_string(const std::string &path, const std::string &value) {
    TrackedDataset dataset;
    dataset.dtype = "object";
    dataset.bytes.assign(value.begin(), value.end());
    datasets_[path] = std::move(dataset);
  }

  void add_bools(const std::string &path,
                 const std::vector<std::size_t> &dimensions,
                 const std::vector<std::uint8_t> &values) {
    TrackedDataset dataset;
    dataset.dtype = "uint8";
    dataset.dimensions = dimensions;
    dataset.bytes = values;
    datasets_[path] = std::move(dataset);
  }

  std::string content_hash(
      const std::string &bundle_file,
      const std::vector<std::string> &path_prefixes = {}) const {
    detail::Sha256 digest;
    digest.update(bundle_file);
    const std::uint8_t separator = 0;
    for (const auto &[path, dataset] : datasets_) {
      bool selected = path_prefixes.empty();
      for (const auto &prefix : path_prefixes)
        selected = selected || path.rfind(prefix, 0) == 0;
      if (!selected)
        continue;
      digest.update(&separator, 1);
      digest.update(path);
      digest.update(&separator, 1);
      digest.update(dataset.dtype);
      digest.update(&separator, 1);
      digest.update(numpy_shape(dataset.dimensions));
      digest.update(&separator, 1);
      if (!dataset.bytes.empty())
        digest.update(dataset.bytes.data(), dataset.bytes.size());
    }
    return "sha256:" + digest.hex_digest();
  }

private:
  std::map<std::string, TrackedDataset> datasets_;
};

class H5File {
public:
  explicit H5File(const std::filesystem::path &path,
                  DatasetHashTracker *tracker = nullptr)
      : handle(path.string(), HighFive::File::Overwrite), tracker(tracker) {}
  H5File(const H5File &) = delete;
  H5File &operator=(const H5File &) = delete;
  HighFive::File handle;
  DatasetHashTracker *tracker;
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
                 const std::vector<std::size_t> &dimensions,
                 const std::vector<T> &values) {
  if (file.tracker)
    file.tracker->add_numeric(path, dimensions, values);
  ensure_groups(file, path);
  auto dataset =
      file.handle.createDataSet<T>(path, HighFive::DataSpace(dimensions));
  if (!values.empty())
    dataset.write_raw(values.data());
}

template <typename T>
void write_scalar(H5File &file, const std::string &path, T value) {
  if (file.tracker)
    file.tracker->add_scalar(path, value);
  ensure_groups(file, path);
  auto dataset =
      file.handle.createDataSet<T>(path, HighFive::DataSpace::From(value));
  dataset.write(value);
}

void write_string(H5File &file, const std::string &path,
                  const std::string &value) {
  if (file.tracker)
    file.tracker->add_string(path, value);
  ensure_groups(file, path);
  auto dataset = file.handle.createDataSet<std::string>(
      path, HighFive::DataSpace::From(value));
  dataset.write(value);
}

void write_bool_array(H5File &file, const std::string &path,
                      const std::vector<std::size_t> &dimensions,
                      const std::vector<std::uint8_t> &values) {
  if (file.tracker)
    file.tracker->add_bools(path, dimensions, values);
  ensure_groups(file, path);
  std::vector<hsize_t> hdf5_dimensions(dimensions.size());
  std::transform(dimensions.begin(), dimensions.end(), hdf5_dimensions.begin(),
                 [](std::size_t dimension) {
                   return static_cast<hsize_t>(dimension);
                 });
  const hid_t space = H5Screate_simple(
      static_cast<int>(hdf5_dimensions.size()), hdf5_dimensions.data(),
      nullptr);
  if (space < 0)
    throw std::runtime_error("failed to create HDF5 boolean dataspace");
  const hid_t datatype = H5Tenum_create(H5T_NATIVE_UCHAR);
  if (datatype < 0) {
    H5Sclose(space);
    throw std::runtime_error("failed to create HDF5 boolean datatype");
  }
  const std::uint8_t false_value = 0;
  const std::uint8_t true_value = 1;
  if (H5Tenum_insert(datatype, "FALSE", &false_value) < 0 ||
      H5Tenum_insert(datatype, "TRUE", &true_value) < 0) {
    H5Tclose(datatype);
    H5Sclose(space);
    throw std::runtime_error("failed to define HDF5 boolean datatype");
  }
  const hid_t dataset = H5Dcreate2(file.handle.getId(), path.c_str(), datatype,
                                   space, H5P_DEFAULT, H5P_DEFAULT,
                                   H5P_DEFAULT);
  if (dataset < 0) {
    H5Tclose(datatype);
    H5Sclose(space);
    throw std::runtime_error("failed to create HDF5 boolean dataset: " + path);
  }
  const bool write_failed =
      !values.empty() &&
      H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
               values.data()) < 0;
  H5Dclose(dataset);
  H5Tclose(datatype);
  H5Sclose(space);
  if (write_failed)
    throw std::runtime_error("failed to write HDF5 boolean dataset: " + path);
}

void write_strings(H5File &file, const std::string &path,
                   const std::vector<std::string> &values) {
  if (file.tracker)
    file.tracker->add_strings(path, values);
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

void finalize_topology(H5File &file, std::size_t atom_count,
                       const std::string &topology_hash,
                       const std::string &atom_order_hash,
                       const std::string &forcefield_hash,
                       const std::string &identity_uuid) {
  const std::string version = kInputSchemaVersion;
  write_string(file, "/schema/name", "sponge.topology.h5");
  write_string(file, "/schema/version", version);
  write_string(file, "/parameters/sponge/schema/name", "sponge.topology.h5");
  write_string(file, "/parameters/sponge/schema/version", version);
  write_string(file, "/identity/uuid", identity_uuid);
  write_string(file, "/topology/atom_order_hash", atom_order_hash);
  write_string(file, "/topology/topology_hash", topology_hash);
  write_string(file, "/topology/forcefield_hash", forcefield_hash);
  write_scalar<std::int64_t>(file, "/topology/atom_count", atom_count);
}

std::pair<float, float> lj_ab(
    const Molecule &molecule, const std::string &lhs,
    const std::string &rhs) {
  const auto a = resolve_molecule_lj_parameter(molecule, lhs);
  const auto b = resolve_molecule_lj_parameter(molecule, rhs);
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

// The legacy serializer interface supplies strings; parse them in memory and
// write typed datasets only. No compatibility files are produced.
struct ListedModule {
  std::string name;
  std::map<std::string, std::string> fields;
};

std::string trim_bundle_text(const std::string &text) {
  const auto begin = text.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos)
    return {};
  return text.substr(begin, text.find_last_not_of(" \t\r\n") - begin + 1);
}

void append_listed_modules(std::vector<ListedModule> &modules,
                           const std::string &text) {
  std::istringstream input(text);
  std::string line, key;
  ListedModule module;
  bool active = false;
  while (std::getline(input, line)) {
    line = trim_bundle_text(line);
    if (line.empty())
      continue;
    if (line.rfind("[[[", 0) == 0 && line.size() >= 6 &&
        line.substr(line.size() - 3) == "]]]") {
      if (active)
        throw std::invalid_argument("listed force missing [[ end ]]");
      module = {trim_bundle_text(line.substr(3, line.size() - 6)), {}};
      if (module.name.empty() ||
          module.name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKL"
                                        "MNOPQRSTUVWXYZ0123456789_") !=
              std::string::npos)
        throw std::invalid_argument("invalid listed force name: " +
                                    module.name);
      active = true;
      key.clear();
    } else if (line.rfind("[[", 0) == 0 && line.size() >= 4 &&
               line.substr(line.size() - 2) == "]]") {
      if (!active)
        throw std::invalid_argument("listed force key outside section");
      key = trim_bundle_text(line.substr(2, line.size() - 4));
      if (key == "end") {
        if (!module.fields.count("potential") ||
            !module.fields.count("parameters") ||
            module.fields.at("potential").empty())
          throw std::invalid_argument(
              "listed force requires potential and parameters: " + module.name);
        auto previous =
            std::find_if(modules.begin(), modules.end(), [&](const auto &item) {
              return item.name == module.name;
            });
        if (previous == modules.end())
          modules.push_back(module);
        else if (previous->fields != module.fields)
          throw std::invalid_argument("conflicting listed force definition: " +
                                      module.name);
        active = false;
      } else {
        if (key != "potential" && key != "parameters" &&
            key != "connected_atoms" && key != "constrain_distance")
          throw std::invalid_argument("unsupported listed force key: " + key);
        if (!module.fields.emplace(key, "").second)
          throw std::invalid_argument("duplicate listed force key: " + key);
      }
    } else {
      if (!active || key.empty())
        throw std::invalid_argument("invalid listed force definition");
      auto &value = module.fields.at(key);
      if (!value.empty())
        value += "\n";
      value += line;
    }
  }
  if (active)
    throw std::invalid_argument("listed force missing [[ end ]]");
}

void write_listed_modules(
    H5File &file, const Molecule &molecule,
    const std::unordered_map<std::string, std::string> &payloads) {
  std::vector<ListedModule> modules;
  auto data = payloads;
  if (!molecule.ryckaert_bellemans.empty()) {
    std::vector<std::pair<std::array<std::int32_t, 4>, std::array<double, 6>>>
        rows;
    for (const auto &term : molecule.ryckaert_bellemans) {
      std::array<double, 6> c{term.c0, term.c1, term.c2,
                              term.c3, term.c4, term.c5};
      if (std::none_of(c.begin(), c.end(), [](double x) { return x != 0; }))
        continue;
      std::array<std::int32_t, 4> atoms{static_cast<std::int32_t>(term.atom0),
                                        static_cast<std::int32_t>(term.atom1),
                                        static_cast<std::int32_t>(term.atom2),
                                        static_cast<std::int32_t>(term.atom3)};
      if (atoms.front() > atoms.back())
        std::reverse(atoms.begin(), atoms.end());
      rows.emplace_back(atoms, c);
    }
    if (!rows.empty()) {
      std::sort(rows.begin(), rows.end());
      modules.push_back({"Ryckaert_Bellemans",
                         {{"parameters", kRyckaertBellemansParameterText},
                          {"potential", kRyckaertBellemansPotential}}});
      std::ostringstream out;
      out << rows.size() << '\n'
          << std::setprecision(std::numeric_limits<double>::max_digits10);
      for (const auto &row : rows) {
        for (auto atom : row.first)
          out << atom << ' ';
        for (auto c : row.second)
          out << c << ' ';
        out << '\n';
      }
      if (!data.emplace("Ryckaert_Bellemans", out.str()).second)
        throw std::invalid_argument("duplicate Ryckaert_Bellemans payload");
    }
  }
  for (const auto &definition : molecule.listed_force_definitions)
    append_listed_modules(modules, definition);
  if (data.count("listed_forces")) {
    append_listed_modules(modules, data.at("listed_forces"));
    data.erase("listed_forces");
  }
  for (const auto &entry : data)
    if (std::none_of(modules.begin(), modules.end(), [&](const auto &module) {
          return module.name == entry.first;
        }))
      throw std::invalid_argument(
          "bundle export does not support active compatibility serializer "
          "without listed force definition: " +
          entry.first);
  if (modules.empty())
    return;
  const std::string root = "/forcefield/custom_force/listed";
  std::vector<std::string> names, potentials, texts, types, parameters,
      connected, distances;
  std::vector<std::int64_t> offsets{0};
  std::vector<std::int32_t> atoms, ints;
  for (const auto &module : modules) {
    const auto &name = module.name;
    std::vector<std::string> local_types, local_names;
    std::vector<std::int32_t> local_atoms;
    std::vector<std::uint8_t> local_ints;
    std::istringstream declarations(module.fields.at("parameters"));
    std::string declaration;
    while (std::getline(declarations, declaration, ',')) {
      std::istringstream parser(declaration);
      std::string type, parameter, extra;
      if (!(parser >> type >> parameter) || (parser >> extra) ||
          (type != "int" && type != "float"))
        throw std::invalid_argument(
            "invalid listed force parameter declaration: " + name);
      if (std::find(local_names.begin(), local_names.end(), parameter) !=
          local_names.end())
        throw std::invalid_argument("duplicate listed force parameter: " +
                                    parameter);
      local_types.push_back(type);
      local_names.push_back(parameter);
      local_ints.push_back(type == "int");
      local_atoms.push_back(type == "int" && parameter.rfind("atom_", 0) == 0 &&
                            parameter.size() == 6);
    }
    const auto atom_parameters =
        std::count(local_atoms.begin(), local_atoms.end(), 1);
    if (atom_parameters < 1 || atom_parameters > 6)
      throw std::invalid_argument(
          "listed force must declare 1 to 6 atom parameters: " + name);
    const auto field = [&](const std::string &key) {
      auto it = module.fields.find(key);
      return it == module.fields.end() ? std::string{} : it->second;
    };
    const auto connection = field("connected_atoms");
    if (!connection.empty() &&
        (connection.size() != 2 || connection[0] == connection[1]))
      throw std::invalid_argument("invalid listed force connected_atoms: " +
                                  name);
    if (!data.count(name))
      throw std::invalid_argument("missing listed force data: " + name);
    std::istringstream input(data.at(name));
    std::int64_t count;
    if (!(input >> count) || count < 0)
      throw std::invalid_argument("invalid listed force item count: " + name);
    std::vector<float> values, floats;
    std::vector<std::int32_t> integers;
    for (std::int64_t i = 0; i < count; ++i) {
      for (std::size_t j = 0; j < local_names.size(); ++j) {
        std::string token;
        if (!(input >> token))
          throw std::invalid_argument("missing listed force parameter data: " +
                                      name);
        std::size_t consumed = 0;
        if (local_ints[j]) {
          const auto value = std::stoll(token, &consumed);
          if (consumed != token.size() ||
              value < std::numeric_limits<std::int32_t>::min() ||
              value > std::numeric_limits<std::int32_t>::max() ||
              (local_atoms[j] &&
               (value < 0 ||
                static_cast<std::size_t>(value) >= molecule.atoms.size())))
            throw std::invalid_argument(
                "invalid listed force integer/atom index: " + name);
          values.push_back(static_cast<float>(value));
          integers.push_back(static_cast<std::int32_t>(value));
          floats.push_back(std::numeric_limits<float>::quiet_NaN());
        } else {
          const auto value = std::stof(token, &consumed);
          if (consumed != token.size() || !std::isfinite(value))
            throw std::invalid_argument("invalid listed force float: " + name);
          values.push_back(value);
          integers.push_back(0);
          floats.push_back(value);
        }
      }
    }
    std::string extra;
    if (input >> extra)
      throw std::invalid_argument("extra listed force parameter data: " + name);
    const auto path = root + "/data/" + name;
    write_string(file, path + "/name", name);
    write_scalar<std::int64_t>(file, path + "/item_count", count);
    write_strings(file, path + "/parameter/name", local_names);
    write_strings(file, path + "/parameter/type", local_types);
    write_bool_array(file, path + "/parameter/is_int", {local_ints.size()},
                     local_ints);
    const std::vector<std::size_t> shape{static_cast<std::size_t>(count),
                                         local_names.size()};
    write_array(file, path + "/parameter/value", shape, values);
    write_array(file, path + "/parameter/int_value", shape, integers);
    write_array(file, path + "/parameter/float_value", shape, floats);
    names.push_back(name);
    potentials.push_back(module.fields.at("potential"));
    texts.push_back(module.fields.at("parameters"));
    connected.push_back(connection);
    distances.push_back(field("constrain_distance"));
    types.insert(types.end(), local_types.begin(), local_types.end());
    parameters.insert(parameters.end(), local_names.begin(), local_names.end());
    atoms.insert(atoms.end(), local_atoms.begin(), local_atoms.end());
    ints.insert(ints.end(), local_ints.begin(), local_ints.end());
    offsets.push_back(parameters.size());
  }
  write_strings(file, root + "/name", names);
  write_strings(file, root + "/potential", potentials);
  write_strings(file, root + "/parameters/text", texts);
  write_strings(file, root + "/parameters/type", types);
  write_strings(file, root + "/parameters/name", parameters);
  write_array(file, root + "/parameters/offset", {offsets.size()}, offsets);
  write_array(file, root + "/parameters/is_atom", {atoms.size()}, atoms);
  write_array(file, root + "/parameters/is_int", {ints.size()}, ints);
  write_strings(file, root + "/connected_atoms", connected);
  write_strings(file, root + "/constrain_distance", distances);
  write_scalar<std::int64_t>(file, root + "/count", modules.size());
}

template <typename Parameter, typename Pair, typename Triple>
void write_manybody_table(
    H5File &file, const Molecule &molecule, const std::string &name,
    const std::string Atom::*atom_field,
    const std::unordered_map<std::string, Parameter> &parameters,
    Pair pair_values, Triple triple_values) {
  if (parameters.empty())
    return;
  std::vector<std::string> types;
  std::vector<std::int32_t> atom_types, pairs, triples;
  std::vector<float> pair_parameters, triple_parameters;
  for (const auto &atom : molecule.atoms) {
    const auto &type = atom.*atom_field;
    if (type.empty())
      throw std::invalid_argument("missing " + name +
                                  " type for atom: " + atom.name);
    auto it = std::find(types.begin(), types.end(), type);
    if (it == types.end()) {
      types.push_back(type);
      it = types.end() - 1;
    }
    atom_types.push_back(static_cast<std::int32_t>(it - types.begin()));
  }
  const auto get = [&](const std::string &key) -> const Parameter & {
    auto it = parameters.find(key);
    if (it == parameters.end())
      throw std::invalid_argument("missing " + name + " parameter: " + key);
    return it->second;
  };
  std::size_t pair_width = 0, triple_width = 0;
  for (std::size_t i = 0; i < types.size(); ++i) {
    for (std::size_t j = 0; j < types.size(); ++j) {
      pairs.insert(pairs.end(), {static_cast<std::int32_t>(i),
                                 static_cast<std::int32_t>(j)});
      const auto pair = pair_values(get(types[i] + "-" + types[j]));
      pair_width = pair.size();
      pair_parameters.insert(pair_parameters.end(), pair.begin(), pair.end());
      for (std::size_t k = 0; k < types.size(); ++k) {
        triples.insert(triples.end(), {static_cast<std::int32_t>(i),
                                       static_cast<std::int32_t>(j),
                                       static_cast<std::int32_t>(k)});
        const auto triple =
            triple_values(get(types[i] + "-" + types[j] + "-" + types[k]));
        triple_width = triple.size();
        triple_parameters.insert(triple_parameters.end(), triple.begin(),
                                 triple.end());
      }
    }
  }
  for (const auto value : pair_parameters)
    if (!std::isfinite(value))
      throw std::invalid_argument("nonfinite " + name + " parameter");
  for (const auto value : triple_parameters)
    if (!std::isfinite(value))
      throw std::invalid_argument("nonfinite " + name + " parameter");
  const auto root = "/manybody/" + name;
  const auto n = types.size();
  write_scalar<std::int32_t>(file, root + "/atom_type_count", n);
  write_array(file, root + "/atom_type", {atom_types.size()}, atom_types);
  write_array(file, root + "/pair/type", {n * n, 2}, pairs);
  write_array(file, root + "/pair/parameters", {n * n, pair_width},
              pair_parameters);
  write_array(file, root + "/triple/type", {n * n * n, 3}, triples);
  write_array(file, root + "/triple/parameters", {n * n * n, triple_width},
              triple_parameters);
}

void write_native_topology(
    H5File &file, const Molecule &molecule,
    const std::unordered_map<std::string, std::string> &listed_payloads) {
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
  const bool materialize_nb14 =
      molecule.write_lj_soft_core || !molecule.nb14_extras.empty();
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
    if (materialize_nb14) {
      const auto lhs =
          resolve_molecule_lj_type(molecule, molecule.atoms[item.atom1].type);
      const auto rhs =
          resolve_molecule_lj_type(molecule, molecule.atoms[item.atom2].type);
      const auto [a, b] = lj_ab(molecule, lhs, rhs);
      nb14_params.insert(nb14_params.end(),
                         {static_cast<float>(item.k_lj * a * 12.0),
                          static_cast<float>(item.k_lj * b * 6.0),
                          static_cast<float>(item.k_ee)});
    } else {
      nb14_params.insert(nb14_params.end(), {static_cast<float>(item.k_lj),
                                             static_cast<float>(item.k_ee)});
    }
  }
  for (const auto &item : molecule.nb14_extras) {
    nb14_atoms.insert(nb14_atoms.end(),
                      {static_cast<std::int32_t>(item.atom1),
                       static_cast<std::int32_t>(item.atom2)});
    nb14_params.insert(nb14_params.end(),
                       {static_cast<float>(item.a * 12.0),
                        static_cast<float>(item.b * 6.0),
                        static_cast<float>(item.kee)});
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
  const std::size_t nb14_count =
      topology.nb14s.size() + molecule.nb14_extras.size();
  write_array(file, "/forcefield/nb14/atoms", {nb14_count, 2}, nb14_atoms);
  write_array(file, "/forcefield/nb14/params",
              {nb14_count, materialize_nb14 ? 3UL : 2UL},
              nb14_params);
  write_scalar<std::int64_t>(file, "/forcefield/nb14/count",
                             nb14_count);

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
    write_scalar<std::int64_t>(file, "/forcefield/urey_bradley/count",
                               molecule.urey_bradleys.size());
  }

  if (!molecule.soft_bonds.empty()) {
    struct SoftBondRow {
      std::int32_t atom1;
      std::int32_t atom2;
      float k;
      float r0;
      std::int32_t from_a_or_b;
    };
    std::vector<SoftBondRow> rows;
    for (const auto &bond : molecule.soft_bonds) {
      if (bond.k == 0.0)
        continue;
      rows.push_back(
          {static_cast<std::int32_t>(std::min(bond.atom1, bond.atom2)),
           static_cast<std::int32_t>(std::max(bond.atom1, bond.atom2)),
           static_cast<float>(bond.k), static_cast<float>(bond.b),
           static_cast<std::int32_t>(bond.from_a_or_b)});
    }
    std::sort(rows.begin(), rows.end(), [](const auto &lhs, const auto &rhs) {
      return std::tie(lhs.atom1, lhs.atom2) < std::tie(rhs.atom1, rhs.atom2);
    });
    if (!rows.empty()) {
      std::vector<std::int32_t> atoms, from_a_or_b;
      std::vector<float> k, r0;
      for (const auto &row : rows) {
        atoms.insert(atoms.end(), {row.atom1, row.atom2});
        k.push_back(row.k);
        r0.push_back(row.r0);
        from_a_or_b.push_back(row.from_a_or_b);
      }
      write_array(file, "/forcefield/bond_soft/atoms", {rows.size(), 2},
                  atoms);
      write_array(file, "/forcefield/bond_soft/k", {rows.size()}, k);
      write_array(file, "/forcefield/bond_soft/r0", {rows.size()}, r0);
      write_array(file, "/forcefield/bond_soft/from_a_or_b", {rows.size()},
                  from_a_or_b);
      write_scalar<std::int64_t>(file, "/forcefield/bond_soft/count",
                                 rows.size());
    }
  }

  write_listed_modules(file, molecule, listed_payloads);
  write_manybody_table(file, molecule, "sw", &Atom::sw_type, molecule.sw_parameters,
      [](const StillingerWeberParameter &p) { return std::array<double, 8>{p.a_big, p.b_big, p.epsilon, p.p, p.q, p.a, p.gamma, p.sigma}; },
      [](const StillingerWeberParameter &p) { return std::array<double, 3>{p.lambda, p.epsilon, p.b}; });
  write_manybody_table(file, molecule, "edip", &Atom::edip_type, molecule.edip_parameters,
      [](const EDIPParameter &p) { return std::array<double, 8>{p.alpha, p.c, p.a, p.a_big, p.b_big, p.rho, p.beta, p.sigma}; },
      [](const EDIPParameter &p) { return std::array<double, 9>{p.eta, p.gamma, p.lambda, p.q0, p.mu, p.u1, p.u2, p.u3, p.u4}; });

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
      const auto name = resolve_molecule_lj_type(molecule, source_type);
      if (!index.count(name)) {
        index[name] = static_cast<std::int32_t>(names.size());
        names.push_back(name);
      }
      atom_type.push_back(index.at(name));
    }
    std::vector<float> pair_a, pair_b;
    for (std::size_t i = 0; i < names.size(); ++i)
      for (std::size_t j = 0; j <= i; ++j) {
        const auto [a, b] = lj_ab(molecule, names[i], names[j]);
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
}

void write_protocol(const std::filesystem::path &path,
                    const std::string &topology_hash,
                    const std::string &protocol_hash,
                    const std::string &identity_uuid) {
  H5File file(path);
  const std::string version = kInputSchemaVersion;
  write_string(file, "/schema/name", "sponge.protocol.h5");
  write_string(file, "/schema/version", version);
  write_string(file, "/parameters/sponge/schema/name", "sponge.protocol.h5");
  write_string(file, "/parameters/sponge/schema/version", version);
  write_string(file, "/identity/uuid", identity_uuid);
  write_string(file, "/protocol/topology_compatibility/topology_hash",
               topology_hash);
  write_string(file, "/identity/content_hash", protocol_hash);
  write_scalar<std::int64_t>(file, "/protocol/cv_count", 0);
  write_scalar<std::int64_t>(file, "/protocol/restraint_count", 0);
}

std::vector<float> restart_box_edges(const Molecule &molecule) {
  const auto box = sponge_coordinate_box_for_export(molecule);
  if (!std::all_of(box.begin(), box.end(),
                   [](double value) { return std::isfinite(value); })) {
    throw std::invalid_argument("box values must be finite");
  }
  const double alpha = box[3] * 3.14159265358979323846 / 180.0;
  const double beta = box[4] * 3.14159265358979323846 / 180.0;
  const double gamma = box[5] * 3.14159265358979323846 / 180.0;
  const double sin_gamma = std::sin(gamma);
  if (std::abs(sin_gamma) < 1e-7) {
    throw std::invalid_argument("box gamma angle produces a singular cell");
  }
  const double cos_alpha = std::cos(alpha);
  const double cos_beta = std::cos(beta);
  const double cos_gamma = std::cos(gamma);
  const double edge_y_x = box[1] * cos_gamma;
  const double edge_y_y = box[1] * sin_gamma;
  const double edge_z_x = box[2] * cos_beta;
  const double edge_z_y =
      box[2] * (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
  const double z2 =
      box[2] * box[2] - edge_z_x * edge_z_x - edge_z_y * edge_z_y;
  if (z2 < -1e-5) {
    throw std::invalid_argument("box angles produce an invalid cell");
  }
  return {
      static_cast<float>(box[0]),
      0,
      0,
      static_cast<float>(edge_y_x),
      static_cast<float>(edge_y_y),
      0,
      static_cast<float>(edge_z_x),
      static_cast<float>(edge_z_y),
      static_cast<float>(std::sqrt(std::max(0.0, z2))),
  };
}

void write_restart(const Molecule &molecule,
                   const std::filesystem::path &path,
                   const std::string &identity_uuid,
                   const std::string &topology_hash,
                   const std::string &atom_order_hash,
                   const std::string &protocol_hash) {
  const auto edges = restart_box_edges(molecule);
  DatasetHashTracker state_tracker;
  H5File file(path, &state_tracker);
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
  write_string(file, "/parameters/sponge/schema/version", kInputSchemaVersion);
  write_string(file, "/schema/name", "sponge.restart.h5");
  write_string(file, "/schema/version", kInputSchemaVersion);
  write_string(file, "/identity/uuid", identity_uuid);
  write_string(file, "/run/topology_hash", topology_hash);
  write_string(file, "/run/atom_order_hash", atom_order_hash);
  write_string(file, "/run/producer_protocol_hash", protocol_hash);
  write_string(
      file, "/run/state_hash",
      state_tracker.content_hash(
          "restart.spgr.h5", {"/particles/", "/parameters/restart/"}));
  write_string(file, "/parameters/sponge/output/status", "finalized");
  write_array<std::int64_t>(file, "/parameters/sponge/output/frame_count", {1},
                            {1});
  write_array<std::int64_t>(
      file, "/parameters/sponge/output/last_complete_step", {1}, {0});
  write_array<double>(file, "/parameters/sponge/output/last_complete_time", {1},
                      {0.0});
  write_array<std::int64_t>(file, "/run/current_step", {1}, {0});
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
                         const std::filesystem::path &dirname,
                         const std::unordered_map<std::string, std::string> &listed_payloads) {
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
  if (molecule.write_min_bonded_parameters) {
    throw std::invalid_argument(
        "bundle export does not support minimum-bonded parameters "
        "(fake_mass, fake_LJ, fake_charge)");
  }
  const std::filesystem::path relative(prefix.empty() ? molecule.name : prefix);
  if (relative.is_absolute())
    throw std::invalid_argument("bundle prefix must be relative");
  std::filesystem::create_directories(dirname);
  const auto base =
      std::filesystem::weakly_canonical(dirname / relative).lexically_normal();
  const auto root = std::filesystem::weakly_canonical(dirname);
  const auto relative_to_root = base.lexically_relative(root);
  if (base == root || relative_to_root.empty() || relative_to_root == "." ||
      *relative_to_root.begin() == "..") {
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

  DatasetHashTracker hash_tracker;
  const std::string identity_uuid = generate_uuid_v4();
  std::string topology_hash;
  std::string atom_order_hash;
  std::string forcefield_hash;
  {
    H5File topology(files.temporary[0], &hash_tracker);
    write_native_topology(topology, molecule, listed_payloads);
    topology_hash = hash_tracker.content_hash("topology.spgt.h5");
    atom_order_hash = hash_tracker.content_hash(
        "topology.spgt.h5", {"/atoms/", "/residues/"});
    forcefield_hash = hash_tracker.content_hash(
        "topology.spgt.h5", {"/forcefield/", "/manybody/", "/qc/"});
    finalize_topology(topology, molecule.atoms.size(), topology_hash,
                      atom_order_hash, forcefield_hash, identity_uuid);
  }
  const DatasetHashTracker empty_protocol_tracker;
  const std::string protocol_hash =
      empty_protocol_tracker.content_hash("protocol.spgp.h5");
  write_protocol(files.temporary[1], topology_hash,
                 protocol_hash, identity_uuid);
  write_restart(molecule, files.temporary[2], identity_uuid, topology_hash,
                atom_order_hash, protocol_hash);
  files.commit();
  return {{"topology", topology_path},
          {"protocol", protocol_path},
          {"restart", restart_path}};
}

} // namespace xpongecpp
