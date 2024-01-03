#include <cstdlib>
#include <exception>
#include <string>

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
    CLI::App app{fmt::format("{} version {}", stmesh::cmake::project_name, stmesh::cmake::project_version)};

    bool show_version = false;
    app.add_flag("--version", show_version, "Show version information");
    CLI11_PARSE(app, argc, argv);

    if (show_version) {
      fmt::print("{}\n", stmesh::cmake::project_version);
      return EXIT_SUCCESS;
    }

    constexpr int five = 5;
    test();
    test_itk();
    test_vtk();
    return factorial(five);
  } catch (const std::exception &e) {
    spdlog::error("Unhandled exception in main: {}", e.what());
  }
}
