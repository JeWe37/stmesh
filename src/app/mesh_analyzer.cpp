#include <cstddef>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <optional>
#include <stmesh/triangulation.hpp>
#include <string>

#include <CLI/App.hpp>
// NOLINTNEXTLINE(misc-include-cleaner)
#include <CLI/CLI.hpp>
#include <CLI/Option.hpp>
#include <CLI/Validators.hpp>
#include <fmt/core.h>
#include <spdlog/cfg/env.h>
#include <spdlog/spdlog.h>

#include <stmesh/stmesh.hpp>

#include <internal_use_only/config.hpp>

int main(int argc, const char **argv) try {
  spdlog::cfg::load_env_levels();
  // NOLINTNEXTLINE(misc-const-correctness)
  CLI::App app{fmt::format("{} version {}", stmesh::cmake::project_name, stmesh::cmake::project_version)};

  // NOLINTBEGIN(misc-const-correctness)
  bool show_version{};
  app.add_flag("--version", show_version, "Show version information");

  std::optional<std::filesystem::path> stats_output_file;
  app.add_option("--statistics-output", stats_output_file, "Write statistics to a file");

  std::optional<std::filesystem::path> vtk_output_dir;
  CLI::Option *vtk_output_dir_option =
      app.add_option("--vtk-output-dir", vtk_output_dir, "Directory to write vtk files to");

  std::string vtk_output_name_format;
  app.add_option("--vtk-output-name-format", vtk_output_name_format, "Format string for vtk output files")
      ->default_val("mesh_{}.vtu")
      ->needs(vtk_output_dir_option);

  // NOLINTBEGIN(*-magic-numbers,cppcoreguidelines-init-variables)
  stmesh::FLOAT_T vtk_output_dt;
  app.add_option("--vtk-output-dt", vtk_output_dt, "Time step for vtk output")
      ->default_val(stmesh::FLOAT_T(0.5))
      ->needs(vtk_output_dir_option);

  std::string vtk_out_coord_format;
  app.add_option("--vtk-out-coord-format", vtk_out_coord_format, "Format string for the coordinates file name")
      ->default_val("")
      ->needs(vtk_output_dir_option);

  std::string vtk_output_vtp_format;
  app.add_option("--vtk-output-vtp-format", vtk_output_vtp_format, "Format string for the vtp output files")
      ->default_val("")
      ->needs(vtk_output_dir_option);

  size_t vtk_output_blocks;
  app.add_option("--vtk-output-blocks", vtk_output_blocks, "Number of blocks for vtk output")
      ->default_val(1)
      ->needs(vtk_output_dir_option);

  stmesh::FLOAT_T output_scale;
  app.add_option("--output-scale", output_scale, "Scale factor for output")->default_val(stmesh::FLOAT_T(1.0));
  // NOLINTEND(*-magic-numbers,cppcoreguidelines-init-variables)

  std::string mixd_file_name;
  app.add_option("mixd_file_name", mixd_file_name, "The name of the MIXD file to read")->required();

  // NOLINTEND(misc-const-correctness)
  CLI11_PARSE(app, argc, argv);

  if (show_version) {
    fmt::print("{}\n", stmesh::cmake::project_version);
    return EXIT_SUCCESS;
  }

  stmesh::TriangulationFromMixd triangulation_from_mixd(mixd_file_name);

  if (vtk_output_dir) {
    spdlog::info("Writing vtk files to {}...", vtk_output_dir->string());
    const stmesh::FLOAT_T min_time = triangulation_from_mixd.boundingBox().min()[3];
    stmesh::writeVTU(*vtk_output_dir, vtk_output_name_format, vtk_output_dt, triangulation_from_mixd, output_scale,
                     min_time, vtk_output_vtp_format, triangulation_from_mixd, vtk_out_coord_format, vtk_output_blocks);
  }
  if (stats_output_file) {
    spdlog::info("Writing statistics to {}...", stats_output_file->string());
    stmesh::writeStatistics(*stats_output_file, triangulation_from_mixd);
  }

  spdlog::info("Done!");
  return EXIT_SUCCESS;
} catch (const std::exception &e) {
  spdlog::error("Unhandled exception in main: {}", e.what());
}