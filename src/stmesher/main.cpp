#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iterator>
#include <optional>
#include <string>

#include <CLI/App.hpp>
// NOLINTNEXTLINE(misc-include-cleaner)
#include <CLI/CLI.hpp>
#include <fmt/core.h>
#include <spdlog/cfg/env.h>
#include <spdlog/spdlog.h>

#include <stmesh/stmesh.hpp>

// This file will be generated automatically when cur_you run the CMake
// configuration step. It creates a namespace called `stmesh`. You can modify
// the source template at `configured_files/config.hpp.in`.
#include <internal_use_only/config.hpp>

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, const char **argv) {
  spdlog::cfg::load_env_levels();
  try {
    // NOLINTNEXTLINE(misc-const-correctness)
    CLI::App app{fmt::format("{} version {}", stmesh::cmake::project_name, stmesh::cmake::project_version)};

    // NOLINTBEGIN(misc-const-correctness)
    bool show_version = false;
    app.add_flag("--version", show_version, "Show version information");

    std::optional<std::filesystem::path> stats_output_file;
    app.add_flag("--statistics-output", stats_output_file, "Write statistics to a file");

    std::optional<std::filesystem::path> vtk_output_dir;
    app.add_option("--vtk-output-dir", vtk_output_dir, "Directory to write vtk files to");

    std::string vtk_output_name_format = "mesh_{}.vtu";
    app.add_option("--vtk-output-name-format", vtk_output_name_format, "Format string for vtk output files");

    // NOLINTBEGIN(*-magic-numbers)
    auto vtk_output_dt = stmesh::FLOAT_T(0.5);
    app.add_option("--vtk-output-dt", vtk_output_dt, "Time step for vtk output");

    auto rho_bar = stmesh::FLOAT_T(20.0);
    app.add_option("--rho-bar", rho_bar, "Rho bar for meshing algorithm");

    auto tau_bar = stmesh::FLOAT_T(0.0013);
    app.add_option("--tau-bar", tau_bar, "Tau bar for meshing algorithm");

    auto zeta = stmesh::FLOAT_T(0.5);
    app.add_option("--zeta", zeta, "Zeta for meshing algorithm");

    auto b = stmesh::FLOAT_T(5.0);
    app.add_option("--b", b, "b for meshing algorithm");

    auto delta = stmesh::FLOAT_T(5.0);
    app.add_option("--delta", delta, "delta for meshing algorithm");
    // NOLINTEND(*-magic-numbers)

    // NOLINTNEXTLINE(cppcoreguidelines-init-variables)
    std::optional<unsigned> seed;
    app.add_option("--seed", seed, "Seed for random number generation");

    // NOLINTEND(misc-const-correctness)
    CLI11_PARSE(app, argc, argv);

    if (show_version) {
      fmt::print("{}\n", stmesh::cmake::project_version);
      return EXIT_SUCCESS;
    }

    const stmesh::SDFSurfaceAdapter<stmesh::HyperSphere4> sdf_surface_adapter(
        stmesh::FLOAT_T(30.0),
        stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(30.0)});

    // NOLINTNEXTLINE(misc-const-correctness)
    stmesh::MeshingAlgorithm meshing_algorithm(sdf_surface_adapter, rho_bar, tau_bar, zeta, b, delta, seed);
    meshing_algorithm.triangulate();

    spdlog::info("Meshing complete! Number of elements: {}",
                 std::distance(meshing_algorithm.triangulation().begin(), meshing_algorithm.triangulation().end()));

    if (vtk_output_dir) {
      spdlog::info("Writing vtk files to {}...", vtk_output_dir->string());
      stmesh::writeVTU(*vtk_output_dir, vtk_output_name_format, vtk_output_dt, sdf_surface_adapter,
                       meshing_algorithm.triangulation());
    }
    if (stats_output_file) {
      spdlog::info("Writing statistics to {}...", stats_output_file->string());
      stmesh::writeStatistics(*stats_output_file, sdf_surface_adapter, meshing_algorithm.triangulation());
    }
    spdlog::info("Done!");
    return EXIT_SUCCESS;
  } catch (const std::exception &e) {
    spdlog::error("Unhandled exception in main: {}", e.what());
  }
}
