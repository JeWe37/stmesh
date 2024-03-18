#include <cstdlib>
#include <exception>
#include <iterator>

#include <CLI/App.hpp>
// NOLINTNEXTLINE(misc-include-cleaner)
#include <CLI/CLI.hpp>
#include <fmt/core.h>
#include <spdlog/spdlog.h>

#include <stmesh/stmesh.hpp>

// This file will be generated automatically when cur_you run the CMake
// configuration step. It creates a namespace called `stmesh`. You can modify
// the source template at `configured_files/config.hpp.in`.
#include <internal_use_only/config.hpp>

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, const char **argv) {
  try {
    // NOLINTNEXTLINE(misc-const-correctness)
    CLI::App app{fmt::format("{} version {}", stmesh::cmake::project_name, stmesh::cmake::project_version)};

    // NOLINTNEXTLINE(misc-const-correctness)
    bool show_version = false;
    app.add_flag("--version", show_version, "Show version information");
    CLI11_PARSE(app, argc, argv);

    if (show_version) {
      fmt::print("{}\n", stmesh::cmake::project_version);
      return EXIT_SUCCESS;
    }

    const stmesh::SDFSurfaceAdapter<stmesh::HyperSphere4> sdf_surface_adapter(
        stmesh::FLOAT_T(30.0),
        stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0)});
    // NOLINTBEGIN(*-magic-numbers,misc-const-correctness)
    stmesh::MeshingAlgorithm meshing_algorithm(sdf_surface_adapter, stmesh::FLOAT_T(20.0), stmesh::FLOAT_T(1e-2),
                                               stmesh::FLOAT_T(0.5), stmesh::FLOAT_T(5.0), stmesh::FLOAT_T(5.0));
    // NOLINTEND(*-magic-numbers,misc-const-correctness)
    meshing_algorithm.triangulate();
    fmt::print("Meshing complete! Number of elements: {}\n",
               std::distance(meshing_algorithm.triangulation().begin(), meshing_algorithm.triangulation().end()));

    return EXIT_SUCCESS;
  } catch (const std::exception &e) {
    spdlog::error("Unhandled exception in main: {}", e.what());
  }
}
