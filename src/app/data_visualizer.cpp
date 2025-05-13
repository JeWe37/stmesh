#include <cstddef>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <stmesh/problem_types.hpp>
#include <stmesh/vtk_writer.hpp>
#include <string>

#include <CLI/App.hpp>
// NOLINTNEXTLINE(misc-include-cleaner)
#include <CLI/CLI.hpp>
#include <CLI/Option.hpp>
#include <CLI/Validators.hpp>
#include <Eigen/Geometry>
#include <fmt/core.h>
#include <spdlog/cfg/env.h>
#include <spdlog/spdlog.h>

#include <stmesh/stmesh.hpp>

#include "transform_subcommand.hpp"

#include <internal_use_only/config.hpp>

int main(int argc, const char **argv) try {
  spdlog::cfg::load_env_levels();
  // NOLINTNEXTLINE(misc-const-correctness)
  CLI::App app{fmt::format("{} version {}(built from {})", stmesh::cmake::project_name, stmesh::cmake::project_version,
                           stmesh::cmake::git_sha)};

  // NOLINTBEGIN(misc-const-correctness)
  CLI::App *version_group = app.add_option_group("Version")->configurable();
  bool show_version{};
  version_group->add_flag("--version", show_version, "Show version information");
  CLI::App *regular_group = app.add_option_group("Regular")->configurable();
  CLI::TriggerOff(version_group, regular_group);

  // NOLINTBEGIN(*-magic-numbers,cppcoreguidelines-init-variables)
  CLI::App *with_coordinates =
      regular_group->add_subcommand("with_coordinates", "Add to vtk files with percomputed coordinates")
          ->configurable();

  std::string vtk_out_coord_format;
  with_coordinates
      ->add_option("vtk-out-coord-format", vtk_out_coord_format, "Format string for the coordinates file name")
      ->required();

  size_t steps;
  with_coordinates->add_option("steps", steps, "Number of steps of mesh")->required();

  CLI::App *without_coordinates =
      regular_group->add_subcommand("without_coordinates", "Add to vtk files without percomputed coordinates")
          ->configurable();

  stmesh::FLOAT_T vtk_output_dt;
  without_coordinates->add_option("--vtk-output-dt", vtk_output_dt, "Time step for vtk output")
      ->default_val(stmesh::FLOAT_T(0.5));

  Eigen::Transform<stmesh::FLOAT_T, 4, Eigen::AffineCompact> transformation =
      Eigen::Transform<stmesh::FLOAT_T, 4, Eigen::AffineCompact>::Identity();
  auto data = addTransformSubcommand(*without_coordinates, transformation);

  stmesh::FLOAT_T min_time{};
  without_coordinates->add_option(
      "--min-time", min_time,
      "The minimum time of the triangulation for offsetting the vtk coordinates when mapping to space-time");

  stmesh::FLOAT_T output_scale{1.0};
  without_coordinates->add_option("--output-scale", output_scale,
                                  "The scale for the vtk coordinates when mapping to space-time");
  // NOLINTEND(*-magic-numbers,cppcoreguidelines-init-variables)

  regular_group->require_subcommand(1);

  std::string mixd_file_name;
  regular_group->add_option("mixd_file_name", mixd_file_name, "The name of the MIXD file to read")->required();

  std::string data_file_name;
  regular_group->add_option("data_file_name", data_file_name, "The name of the data file to read")->required();

  std::string problem_type_str;
  regular_group->add_option("problem_type", problem_type_str, "The problem type to use for the data file")
      ->required()
      ->check(CLI::IsMember(stmesh::kNameMap));

  std::filesystem::path vtk_output_dir;
  regular_group->add_option("vtk-output-dir", vtk_output_dir, "Directory to write vtk files to")
      ->required()
      ->check(CLI::ExistingDirectory);

  std::string vtk_mesh_name_format;
  regular_group->add_option("--vtk-mesh-name-format", vtk_mesh_name_format, "Format string for vtk mesh files")
      ->default_val("mesh_{}.vtu");

  std::string vtk_output_name_format;
  regular_group->add_option("--vtk-output-name-format", vtk_output_name_format, "Format string for vtk mesh files")
      ->default_val("simulation_{}.vtu");

  std::string write_config;
  regular_group->add_option("--write-config", write_config, "Write the config to a file");

  app.set_config("--config")->required(false);
  app.callback([&]() { spdlog::info("Using config:\n{}", app.config_to_str()); });

  // NOLINTEND(misc-const-correctness)
  CLI11_PARSE(app, argc, argv);

  if (show_version) {
    fmt::print("{} commit {}\n", stmesh::cmake::project_version, stmesh::cmake::git_sha);
    return EXIT_SUCCESS;
  }
  if (!write_config.empty()) {
    auto formatter = app.get_config_formatter();
    std::ofstream config_file(write_config);
    config_file << formatter->to_config(&app, true, true, "");
    return EXIT_SUCCESS;
  }

  const stmesh::TriangulationFromMixdWithData triangulation_from_mixd(
      mixd_file_name, stmesh::kNameMap.at(problem_type_str), data_file_name);

  if (regular_group->got_subcommand(with_coordinates))
    stmesh::addVTUData(vtk_output_dir, vtk_mesh_name_format, vtk_output_name_format, steps, triangulation_from_mixd,
                       vtk_out_coord_format);
  else
    stmesh::addVTUData(vtk_output_dir, vtk_mesh_name_format, vtk_output_name_format, vtk_output_dt,
                       triangulation_from_mixd, transformation, output_scale, min_time);

  spdlog::info("Done!");
  return EXIT_SUCCESS;
} catch (const std::exception &e) {
  spdlog::error("Unhandled exception in main: {}", e.what());
}