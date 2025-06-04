#ifndef STMESH_TRANSFORM_SUBCOMMAND_HPP
#define STMESH_TRANSFORM_SUBCOMMAND_HPP
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include <CLI/App.hpp>
#include <CLI/CLI.hpp>
#include <CLI/Option.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <stmesh/stmesh.hpp>

struct TransformData {
  std::array<stmesh::FLOAT_T, 4> translation{};
  stmesh::FLOAT_T rotate_xy{};
  stmesh::FLOAT_T rotate_xz{};
  stmesh::FLOAT_T rotate_xw{};
  stmesh::FLOAT_T rotate_yz{};
  stmesh::FLOAT_T rotate_yw{};
  stmesh::FLOAT_T rotate_zw{};
  std::array<stmesh::FLOAT_T, 4> scale{};
  std::vector<stmesh::FLOAT_T> custom_matrix;
};

[[nodiscard]] inline TransformData
addTransformSubcommand(CLI::App &app, Eigen::Transform<stmesh::FLOAT_T, 4, Eigen::AffineCompact> &transform_matrix) {
  TransformData data;
  auto *transform = app.add_subcommand("transform", "Apply 4D affine transformation");

  transform_matrix.setIdentity();

  // Translation
  CLI::Option *translation_opt =
      transform->add_option("-t,--translate", data.translation, "4D translation vector (x,y,z,w)")->expected(4);

  // Individual rotation angles (in radians)
  CLI::Option *rotation_xy_opt =
      transform->add_option("--rotate-xy", data.rotate_xy, "Rotation angle in radians for XY plane");

  CLI::Option *rotation_xz_opt =
      transform->add_option("--rotate-xz", data.rotate_xz, "Rotation angle in radians for XZ plane");

  CLI::Option *rotation_xw_opt =
      transform->add_option("--rotate-xw", data.rotate_xw, "Rotation angle in radians for XW plane");

  CLI::Option *rotation_yz_opt =
      transform->add_option("--rotate-yz", data.rotate_yz, "Rotation angle in radians for YZ plane");

  CLI::Option *rotation_yw_opt =
      transform->add_option("--rotate-yw", data.rotate_yw, "Rotation angle in radians for YW plane");

  CLI::Option *rotation_zw_opt =
      transform->add_option("--rotate-zw", data.rotate_zw, "Rotation angle in radians for ZW plane");

  // Scale factors
  std::ranges::fill(data.scale, 1.0);
  CLI::Option *scale_opt = transform->add_option("-s,--scale", data.scale, "Scale factors (x,y,z,w)")->expected(4);

  // Custom matrix elements (row-major)
  transform->add_option("--matrix", data.custom_matrix, "Custom transform matrix (20 elements, row-major)")
      ->expected(4 * 5)
      ->excludes(translation_opt)
      ->excludes(rotation_xy_opt)
      ->excludes(rotation_xz_opt)
      ->excludes(rotation_xw_opt)
      ->excludes(rotation_yz_opt)
      ->excludes(rotation_yw_opt)
      ->excludes(rotation_zw_opt)
      ->excludes(scale_opt);

  transform->callback([&]() {
    // Override with custom matrix if provided
    if (!data.custom_matrix.empty()) {
      transform_matrix.matrix() = Eigen::Map<Eigen::Matrix<stmesh::FLOAT_T, 4, 5>>{data.custom_matrix.data()};
      return;
    }

    // Apply translation
    transform_matrix.translate(
        stmesh::Vector4F(data.translation[0], data.translation[1], data.translation[2], data.translation[3]));

    // Helper function for rotation matrix creation
    auto create_rotation_matrix = [](int plane1, int plane2, stmesh::FLOAT_T angle) {
      Eigen::Matrix<stmesh::FLOAT_T, 4, 4> rot = Eigen::Matrix<stmesh::FLOAT_T, 4, 4>::Identity();
      stmesh::FLOAT_T c = std::cos(angle);
      stmesh::FLOAT_T s = std::sin(angle);
      rot(plane1, plane1) = c;
      rot(plane1, plane2) = -s;
      rot(plane2, plane1) = s;
      rot(plane2, plane2) = c;
      return rot;
    };

    // Apply rotations
    transform_matrix *= create_rotation_matrix(0, 1, data.rotate_xy) * create_rotation_matrix(0, 2, data.rotate_xz) *
                        create_rotation_matrix(0, 3, data.rotate_xw) * create_rotation_matrix(1, 2, data.rotate_yz) *
                        create_rotation_matrix(1, 3, data.rotate_yw) * create_rotation_matrix(2, 3, data.rotate_zw);

    // Apply scaling
    transform_matrix.scale(Eigen::Map<stmesh::Vector4F>(data.scale.data()));
  });

  return data;
}
#endif