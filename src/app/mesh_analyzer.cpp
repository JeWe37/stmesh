#include <cstddef>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <optional>
#include <stmesh/problem_types.hpp>
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

  std::optional<std::filesystem::path> stats_output_file;
  regular_group->add_option("--statistics-output", stats_output_file, "Write statistics to a file");

  std::optional<std::filesystem::path> vtk_output_dir;
  CLI::Option *vtk_output_dir_option =
      regular_group->add_option("--vtk-output-dir", vtk_output_dir, "Directory to write vtk files to")
          ->check(CLI::ExistingDirectory);

  std::string vtk_output_name_format;
  regular_group->add_option("--vtk-output-name-format", vtk_output_name_format, "Format string for vtk output files")
      ->default_val("mesh_{}.vtu")
      ->needs(vtk_output_dir_option);

  // NOLINTBEGIN(*-magic-numbers,cppcoreguidelines-init-variables)
  stmesh::FLOAT_T vtk_output_dt;
  regular_group->add_option("--vtk-output-dt", vtk_output_dt, "Time step for vtk output")
      ->default_val(stmesh::FLOAT_T(0.5))
      ->needs(vtk_output_dir_option);

  std::string vtk_out_coord_format;
  regular_group
      ->add_option("--vtk-out-coord-format", vtk_out_coord_format, "Format string for the coordinates file name")
      ->default_val("")
      ->needs(vtk_output_dir_option);

  std::string vtk_output_vtp_format;
  regular_group->add_option("--vtk-output-vtp-format", vtk_output_vtp_format, "Format string for the vtp output files")
      ->default_val("")
      ->needs(vtk_output_dir_option);

  size_t vtk_output_blocks;
  regular_group->add_option("--vtk-output-blocks", vtk_output_blocks, "Number of blocks for vtk output")
      ->default_val(1)
      ->needs(vtk_output_dir_option);

  stmesh::FLOAT_T output_scale;
  regular_group->add_option("--output-scale", output_scale, "Scale factor for output")
      ->default_val(stmesh::FLOAT_T(1.0));
  // NOLINTEND(*-magic-numbers,cppcoreguidelines-init-variables)

  std::optional<std::string> data_file;
  CLI::Option *data_file_option =
      regular_group->add_option("--data-file", data_file, "The name of the data file to read");

  std::string problem_type_str;
  data_file_option->needs(
      regular_group->add_option("--problem-type", problem_type_str, "The problem type to use for the data file")
          ->check(CLI::IsMember(stmesh::kNameMap))
          ->needs(data_file_option));

  std::string mixd_file_name;
  regular_group->add_option("mixd_file_name", mixd_file_name, "The name of the MIXD file to read")->required();

  Eigen::Transform<stmesh::FLOAT_T, 4, Eigen::AffineCompact> transformation =
      Eigen::Transform<stmesh::FLOAT_T, 4, Eigen::AffineCompact>::Identity();
  auto data = addTransformSubcommand(app, transformation);

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

  auto output = [&](const auto &triangulation_from_mixd) {
    if (vtk_output_dir) {
      spdlog::info("Writing vtk files to {}...", vtk_output_dir->string());
      const stmesh::FLOAT_T min_time = triangulation_from_mixd.boundingBox().min()[3];
      stmesh::writeVTU(*vtk_output_dir, vtk_output_name_format, vtk_output_dt, triangulation_from_mixd, transformation,
                       output_scale, min_time, vtk_output_vtp_format, triangulation_from_mixd, vtk_out_coord_format,
                       vtk_output_blocks);
    }
    if (stats_output_file) {
      spdlog::info("Writing statistics to {}...", stats_output_file->string());
      stmesh::writeStatistics(*stats_output_file, triangulation_from_mixd);
    }
  };

  if (data_file) {
    spdlog::info("Adding data from {}...", *data_file);
    const stmesh::TriangulationFromMixdWithData triangulation_from_mixd(
        mixd_file_name, stmesh::kNameMap.at(problem_type_str), *data_file);
    output(triangulation_from_mixd);
  } else {
    const stmesh::TriangulationFromMixd triangulation_from_mixd(mixd_file_name);
    output(triangulation_from_mixd);
  }

  spdlog::info("Done!");
  return EXIT_SUCCESS;
} catch (const std::exception &e) {
  spdlog::error("Unhandled exception in main: {}", e.what());
}