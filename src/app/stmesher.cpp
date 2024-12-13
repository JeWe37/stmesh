#include <array>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <stmesh/lfs_schemes.hpp>
#include <stmesh/utility.hpp>
#include <string>
#include <vector>

#include <CLI/App.hpp>
// NOLINTNEXTLINE(misc-include-cleaner)
#include <CLI/CLI.hpp>
#include <CLI/Option.hpp>
#include <CLI/Validators.hpp>
#include <exprtk.hpp>
#include <fmt/core.h>
#include <spdlog/cfg/env.h>
#include <spdlog/spdlog.h>

#include <stmesh/stmesh.hpp>

// This file will be generated automatically when cur_you run the CMake
// configuration step. It creates a namespace called `stmesh`. You can modify
// the source template at `configured_files/config.hpp.in`.
#include <internal_use_only/config.hpp>

enum class RadiusSchemes : uint8_t { kConstant, kImage, kLfs, kBoundary };

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

  stmesh::FLOAT_T rho_bar;
  app.add_option("--rho-bar", rho_bar, "Rho bar for meshing algorithm")->default_val(stmesh::FLOAT_T(20.0));

  stmesh::FLOAT_T tau_bar;
  app.add_option("--tau-bar", tau_bar, "Tau bar for meshing algorithm")->default_val(stmesh::FLOAT_T(0.0013));

  stmesh::FLOAT_T zeta;
  app.add_option("--zeta", zeta, "Zeta for meshing algorithm")->default_val(stmesh::FLOAT_T(0.5));

  stmesh::FLOAT_T b;
  app.add_option("--b", b, "b for meshing algorithm")->default_val(stmesh::FLOAT_T(5.0));

  stmesh::FLOAT_T delta;
  app.add_option("--delta", delta, "delta for meshing algorithm")->default_val(stmesh::FLOAT_T(5.0));

  RadiusSchemes radius_scheme{};
  std::map<std::string, RadiusSchemes> radius_scheme_map{{"constant", RadiusSchemes::kConstant},
                                                         {"image", RadiusSchemes::kImage},
                                                         {"lfs", RadiusSchemes::kLfs},
                                                         {"boundary", RadiusSchemes::kBoundary}};
  app.add_option_function<std::string>(
         "--radius-scheme", [&](const std::string &val) { radius_scheme = radius_scheme_map.at(val); },
         "Radius scheme for meshing algorithm")
      ->default_val("constant")
      ->run_callback_for_default();

  std::string radius_scheme_arg;
  app.add_option("--radius-scheme-arg", radius_scheme_arg, "Argument for the radius scheme")
      ->default_val(std::numeric_limits<stmesh::FLOAT_T>::infinity());

  bool disable_rule6{};
  app.add_flag("--disable-rule6", disable_rule6, "Disable picking region");

  std::optional<unsigned> seed;
  app.add_option("--seed", seed, "Seed for random number generation");
  // NOLINTEND(*-magic-numbers,cppcoreguidelines-init-variables)

  std::optional<std::string> edt_file;
  CLI::Option *edt_file_option = app.add_option("--edt-file", edt_file, "Read an EDT file");

  bool constant_lfs = true;
  app.add_flag("!--no-constant-lfs", constant_lfs, "Use a constant local feature size, default true")
      ->needs(edt_file_option);

  stmesh::HypercubeBoundaryManager hypercube_boundary_manager;
  CLI::Option *hypercube_option =
      app.add_option_function<std::vector<std::array<stmesh::FLOAT_T, 2ULL * 4ULL>>>(
             "--hypercube",
             [&hypercube_boundary_manager](const auto &vec) {
               for (const auto &arr : vec) {
                 stmesh::Vector4F min;
                 std::copy_n(arr.begin(), 4, min.data());
                 stmesh::Vector4F max;
                 std::copy_n(arr.begin() + 4, 4, max.data());

                 stmesh::HyperCube4 hypercube(min, max);
                 fmt::print("Added hypercube {} with id {}\n", hypercube,
                            static_cast<int>(hypercube_boundary_manager.addBoundaryRegion(hypercube)));
               }
             },
             "Add a hypercube to the meshing algorithm")
          ->take_all();

  bool use_edt_file_boundary_regions = false;
  app.add_flag("--use-edt-file-boundary-regions", use_edt_file_boundary_regions,
               "Use boundary regions from the EDT file")
      ->needs(edt_file_option)
      ->excludes(hypercube_option);

  std::optional<std::filesystem::path> mixd_output_file;
  app.add_option("--mixd-output", mixd_output_file,
                 "Specify the .minf file to write, other MIXD files will be placed alongside it.");

  // NOLINTEND(misc-const-correctness)
  CLI11_PARSE(app, argc, argv);

  if (show_version) {
    fmt::print("{}\n", stmesh::cmake::project_version);
    return EXIT_SUCCESS;
  }

  const auto mesh = [&]<bool LFS = false>(const auto &surface_adapter, auto &&lfs_scheme, auto &&selected_radius_scheme,
                                          const std::shared_ptr<stmesh::EDTReader<4, LFS>> &edt_reader = nullptr) {
    // NOLINTNEXTLINE(misc-const-correctness)
    stmesh::MeshingAlgorithm meshing_algorithm(
        surface_adapter, rho_bar, tau_bar, zeta, b, std::forward<decltype(lfs_scheme)>(lfs_scheme),
        std::forward<decltype(selected_radius_scheme)>(selected_radius_scheme), seed, disable_rule6);

    spdlog::info("Setup complete. Starting meshing...");

    meshing_algorithm.triangulate();

    spdlog::info("Meshing complete! Number of elements: {}",
                 std::distance(meshing_algorithm.triangulation().begin(), meshing_algorithm.triangulation().end()));

    const stmesh::FLOAT_T min_time =
        (!vtk_out_coord_format.empty() || mixd_output_file) ? meshing_algorithm.minTime() : stmesh::FLOAT_T();

    spdlog::debug("Minimum time in mesh: {}", min_time);

    auto writable_triangulation = meshing_algorithm.triangulation().writableTriangulation(&surface_adapter);
    if (vtk_output_dir) {
      spdlog::info("Writing vtk files to {}...", vtk_output_dir->string());
      if (use_edt_file_boundary_regions)
        stmesh::writeVTU(*vtk_output_dir, vtk_output_name_format, vtk_output_dt, writable_triangulation, output_scale,
                         min_time, vtk_output_vtp_format, *edt_reader, vtk_out_coord_format, vtk_output_blocks);
      else
        stmesh::writeVTU(*vtk_output_dir, vtk_output_name_format, vtk_output_dt, writable_triangulation, output_scale,
                         min_time, vtk_output_vtp_format, hypercube_boundary_manager, vtk_out_coord_format,
                         vtk_output_blocks);
    }
    if (mixd_output_file) {
      spdlog::info("Writing MIXD files to {}...", mixd_output_file->string());
      if (use_edt_file_boundary_regions)
        stmesh::writeMixd(*mixd_output_file, surface_adapter, writable_triangulation, *edt_reader, output_scale,
                          min_time);
      else
        stmesh::writeMixd(*mixd_output_file, surface_adapter, writable_triangulation, hypercube_boundary_manager,
                          output_scale, min_time);
    }
    if (stats_output_file) {
      spdlog::info("Writing statistics to {}...", stats_output_file->string());
      stmesh::writeStatistics(*stats_output_file, writable_triangulation);
    }
  };

  const auto invoke_scheme = [&]<bool LFS = false>(const auto &surface_adapter, auto &&lfs_scheme,
                                                   const std::shared_ptr<stmesh::EDTReader<4, LFS>> &edt_reader =
                                                       nullptr) {
    if (radius_scheme == RadiusSchemes::kConstant) {
      spdlog::debug("Using constant radius scheme");
      mesh(surface_adapter, std::forward<decltype(lfs_scheme)>(lfs_scheme),
           stmesh::radius_schemes::Constant(stmesh::FLOAT_T(std::stod(radius_scheme_arg))), edt_reader);
    } else if (radius_scheme == RadiusSchemes::kImage) {
      spdlog::debug("Using image radius scheme");
      mesh(surface_adapter, std::forward<decltype(lfs_scheme)>(lfs_scheme),
           stmesh::radius_schemes::ImageRadius(radius_scheme_arg), edt_reader);
    } else {
      // NOLINTBEGIN(misc-const-correctness)
      exprtk::symbol_table<stmesh::FLOAT_T> symbol_table;
      exprtk::expression<stmesh::FLOAT_T> expression;
      exprtk::parser<stmesh::FLOAT_T> parser;

      stmesh::FLOAT_T d{};
      // NOLINTEND(misc-const-correctness)
      symbol_table.add_variable("d", d);
      symbol_table.add_constants();

      expression.register_symbol_table(symbol_table);

      if (!parser.compile(radius_scheme_arg, expression)) {
        spdlog::error("Failed to compile expression: {}", parser.error());
        return false;
      }

      auto radius_scheme_lambda = [&](const stmesh::FLOAT_T &val) {
        d = val;
        return expression.value();
      };

      if (radius_scheme == RadiusSchemes::kLfs) {
        if constexpr (LFS) {
          spdlog::debug("Using LFS radius scheme");
          mesh(surface_adapter, std::forward<decltype(lfs_scheme)>(lfs_scheme),
               stmesh::radius_schemes::LFSRadius(edt_reader, radius_scheme_lambda), edt_reader);
        } else
          return false;
      } else if (radius_scheme == RadiusSchemes::kBoundary && edt_reader) {
        spdlog::debug("Using boundary distance radius scheme");
        mesh(surface_adapter, std::forward<decltype(lfs_scheme)>(lfs_scheme),
             stmesh::radius_schemes::BoundaryDistanceRadius(edt_reader, radius_scheme_lambda), edt_reader);
      } else
        return false;
    }
    return true;
  };

  if (edt_file) {
    spdlog::info("Reading EDT file {}...", *edt_file);
    if (constant_lfs) {
      spdlog::debug("Using constant LFS");
      const auto edt_reader = std::make_shared<stmesh::EDTReader<4>>(*edt_file);
      const stmesh::EDTSurfaceAdapter adapter(edt_reader);
      spdlog::info("EDT file read successfully! Bounding box: min=({}), max=({})",
                   edt_reader->boundingBox().min().transpose(), edt_reader->boundingBox().max().transpose());
      if (!invoke_scheme(adapter, stmesh::lfs_schemes::Constant(delta), edt_reader))
        return EXIT_FAILURE;
    } else {
      spdlog::debug("Using binary image approximation LFS");
      const auto edt_reader = std::make_shared<stmesh::EDTReader<4, true>>(*edt_file);
      const stmesh::EDTSurfaceAdapter adapter(edt_reader);
      spdlog::info("EDT file read successfully! Bounding box: min=({}), max=({})",
                   edt_reader->boundingBox().min().transpose(), edt_reader->boundingBox().max().transpose());
      if (!invoke_scheme(adapter, stmesh::lfs_schemes::BinaryImageApproximation(delta, edt_reader), edt_reader))
        return EXIT_FAILURE;
    }
  } else {
    const stmesh::SDFSurfaceAdapter<stmesh::HyperSphere4> sdf_surface_adapter(
        stmesh::FLOAT_T(30.0),
        stmesh::Vector4F{stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(0.0), stmesh::FLOAT_T(30.0)});
    spdlog::debug("Using constant LFS");
    if (!invoke_scheme(sdf_surface_adapter, stmesh::lfs_schemes::Constant(delta)))
      return EXIT_FAILURE;
  }
  spdlog::info("Done!");
  return EXIT_SUCCESS;
} catch (const std::exception &e) {
  spdlog::error("Unhandled exception in main: {}", e.what());
}
