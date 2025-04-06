#ifndef STMESH_PROBLEM_TYPES_HPP
#define STMESH_PROBLEM_TYPES_HPP

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <map>
#include <numeric>
#include <span>
#include <string>
#include <vector>

#include "utility.hpp"

namespace stmesh {
/// A data entry for a problem type
struct DataEntry {
  const char *name; ///< The name of the data entry
  size_t length;    ///< The length of the data entry
};

/// A problem type
/**
 * A problem type is a collection of data entries. Each data entry has a name and a length. The problem type is used to
 * split data into its constituent parts.
 */
struct ProblemType {
  std::vector<DataEntry> data_entries; ///< The data entries of the problem type

  /// Get the number of entries in the problem type
  /**
   * Get the number of entries in the problem type. This matches the expected size of the data vector for the problem.
   *
   * @return The number of entries in the problem type
   */
  [[nodiscard]] constexpr size_t entries() const noexcept {
    return std::accumulate(data_entries.begin(), data_entries.end(), size_t{0},
                           [](size_t acc, const DataEntry &entry) { return acc + entry.length; });
  }

  /// Loop over the data entries
  /**
   * Loop over the data entries. This function loops over the data entries in the problem type, calling the function f
   * with the index of the data entry and a span of the data for that entry.
   *
   * @tparam F The type of the function to call
   * @param data The data to loop over
   * @param f The function to call
   */
  template <typename F>
  constexpr void forEach(const Eigen::Vector<FLOAT_T, Eigen::Dynamic> &data, const F &f) const
  requires std::invocable<F, size_t, std::span<const FLOAT_T>>
  {
    size_t index = 0;
    for (size_t i = 0; i < data_entries.size(); ++i) {
      f(i, std::span<const FLOAT_T>(&data[static_cast<Eigen::Index>(index)], data_entries[i].length));
      index += data_entries[i].length;
    }
  }
};

// NOLINTBEGIN(cert-err58-cpp)
const ProblemType kSolidProblem = {{{"displacement", 3}}}; ///< The problem type for a solid problem
const ProblemType kViscoelasticProblem = {
    {DataEntry{"velocity", 3}, DataEntry{"stress", 6}}};                          ///< The problem type for a
                                                                                  ///< viscoelastic problem
const ProblemType kAdvectionDiffusionProblem = {{DataEntry{"concentration", 1}}}; ///< The problem type for an advection
                                                                                  ///< diffusion problem
const ProblemType kINSProblem = {
    {DataEntry{"velocity", 3}, DataEntry{"pressure", 1}}}; ///< The problem type for an incompressible Navier-Stokes
                                                           ///< problem
const ProblemType kCNSProblem = {{DataEntry{"pressure", 1}, DataEntry{"velocity", 3},
                                  DataEntry{"temperature", 1}}}; ///< The problem type for a compressible Navier-Stokes
                                                                 ///< problem
const ProblemType kRVTCNSProblem = {{DataEntry{"density", 1}, DataEntry{"velocity", 3},
                                     DataEntry{"temperature", 1}}}; ///< The problem type for a density based
                                                                    ///< compressible Navier-Stokes problem

#pragma GCC diagnostic push // see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=55776
#pragma GCC diagnostic ignored "-Wshadow"
enum class ProblemTypeEnum : std::uint8_t {
  kSolidProblem,
  kViscoelasticProblem,
  kAdvectionDiffusionProblem,
  kINSProblem,
  kCNSProblem,
  kRVTCNSProblem
}; ///< The problem type enum
#pragma GCC diagnostic pop

const std::map<ProblemTypeEnum, ProblemType> kProblemTypeMap = {
    {ProblemTypeEnum::kSolidProblem, kSolidProblem},
    {ProblemTypeEnum::kViscoelasticProblem, kViscoelasticProblem},
    {ProblemTypeEnum::kAdvectionDiffusionProblem, kAdvectionDiffusionProblem},
    {ProblemTypeEnum::kINSProblem, kINSProblem},
    {ProblemTypeEnum::kCNSProblem, kCNSProblem},
    {ProblemTypeEnum::kRVTCNSProblem, kRVTCNSProblem}}; ///< The map from problem type enums to problem types

const std::map<std::string, ProblemType> kNameMap = {{"solid", kSolidProblem},
                                                     {"viscoelastic", kViscoelasticProblem},
                                                     {"advection_diffusion", kAdvectionDiffusionProblem},
                                                     {"ins", kINSProblem},
                                                     {"cns", kCNSProblem},
                                                     {"rvtcns", kRVTCNSProblem}}; ///< The map from problem type names
                                                                                  ///< to problem types
// NOLINTEND(cert-err58-cpp)
} // namespace stmesh

#endif