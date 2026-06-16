#pragma once

#include <vector>

#include <Eigen/Dense>

#include "config.hpp"
#include "types.hpp"

namespace weather {

// Energy (RMS) average of a stack of hourly scalar fields:
//   H_eff(i,j) = sqrt( mean_{t valid} H_t(i,j)^2 )
// A cell is NaN in the input snapshot if and only if that hour's value at
// that cell is missing (land/ice). A cell in the output is NaN if it has
// fewer than `min_valid_hours` valid (non-NaN) input timesteps.
// Used for sigwh, wsh (height fields) and was (drag-equivalent wind speed).
Eigen::MatrixXd energy_average(const EnergyStack& stack,
                                int min_valid_hours = config::DEFAULT_MIN_VALID_HOURS);

// Energy-flux-weighted mean period:
//   T_eff(i,j) = sum_{t valid} H_t^2 * T_t  /  sum_{t valid} H_t^2
// A timestep contributes weight H_t^2 to the denominator whenever H_t is
// valid; it additionally contributes to the numerator only when T_t is also
// valid. Output is NaN where the cell has fewer than `min_valid_hours` valid
// height timesteps, or where the weight sum is zero (a calm day).
Eigen::MatrixXd energy_weighted_period(
    const EnergyStack& height_stack, const EnergyStack& period_stack,
    int min_valid_hours = config::DEFAULT_MIN_VALID_HOURS);

// Energy-weighted circular mean of a directional field (meteorological
// convention: clockwise from North):
//   u = sum H_t^2 sin(theta_t),  v = sum H_t^2 cos(theta_t)
//   theta_eff = atan2(u, v) mod 360
// A timestep contributes only when both the height and the angle at that
// timestep are valid. Output is NaN where the cell has fewer than
// `min_valid_hours` such jointly-valid timesteps.
Eigen::MatrixXd energy_weighted_circular_mean(
    const EnergyStack& height_stack, const EnergyStack& angle_deg_stack,
    int min_valid_hours = config::DEFAULT_MIN_VALID_HOURS);

// Non-local (swell) wave energy residual, derived from already daily-averaged
// fields:
//   H_swell(i,j) = sqrt( max( sigwh_daily(i,j)^2 - wsh_daily(i,j)^2, 0 ) )
// NaN propagates if either input is NaN. Floating point noise that would
// otherwise produce a small negative radicand is clamped to zero.
Eigen::MatrixXd swell_residual(const Eigen::MatrixXd& sigwh_daily,
                                const Eigen::MatrixXd& wsh_daily);

// Combines one day's worth of hourly snapshots (in chronological order) into
// the 6 daily-averaged output fields described in README.md Section 3.
DailyAverages daily_energy_average(
    const std::vector<WindSeaResult>& hourly_snapshots,
    int min_valid_hours = config::DEFAULT_MIN_VALID_HOURS);

// Decodes a 16-point compass-rose direction code into degrees (clockwise
// from North), per the GFS-Wave archive's int8-coded direction fields:
//   degrees = (code - origin_code) * 22.5,  wrapped into [0, 360)
// `origin_code` is the code that represents North -- 0 for wsd (codes 0-15)
// and 1 for wad (codes 1-16) in the archive observed in practice. NaN
// propagates unchanged.
Eigen::MatrixXd compass_code_to_degrees(const Eigen::MatrixXd& codes,
                                         double origin_code);

} // namespace weather
