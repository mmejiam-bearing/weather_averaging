#pragma once

#include <array>

namespace weather::config {

// The real GFS-Wave archive splits fields across two grids rather than one
// uniform 721x1440 global grid (see README.md Section 2).
inline constexpr int WAVE_GRID_ROWS = 621; // sigwh, wsh, wsp, wsd, swell_residual
inline constexpr int WIND_GRID_ROWS = 721; // was, wad
inline constexpr int GRID_COLS = 1440;     // shared by both grids

// Latitude mapping: lat(row) = LAT_AT_ROW0 + row * LAT_STEP_DEG.
//
// Both grids store row 0 at the SOUTH end, with latitude increasing as row
// increases (northward) -- this was confirmed two independent ways:
//   - wave grid: correlating sigwh's NaN (land/ice) mask against Natural
//     Earth coastlines (92.6% agreement at a single sharp peak, verified
//     visually -- continents line up with coastlines exactly).
//   - wind grid: `was` carries no NaN/land mask to correlate against, so it
//     was verified by (a) the user visually confirming the same Patagonia/
//     South Africa displacement artifact seen in the wave grid's raw (no
//     map) render, and (b) independently, strong winter-storm signatures in
//     `was` landing at the row a ~45-50N North Pacific/Atlantic cyclone
//     should occupy under this mapping (and an implausible Southern Ocean
//     latitude under the opposite one).
// An earlier version of this pipeline assumed the wind grid followed the
// opposite (row 0 = north, decreasing) textbook convention; that assumption
// was never independently verified and was wrong -- do not reintroduce it.
inline constexpr double WAVE_GRID_LAT_AT_ROW0 = -75.0; // wave grid row 0 (south end)
inline constexpr double WIND_GRID_LAT_AT_ROW0 = -90.0; // wind grid row 0 (south pole)
inline constexpr double GRID_LAT_STEP_DEG = 0.25;       // both grids: lat increases northward with row
inline constexpr double GRID_LON_RESOLUTION_DEG = 0.25; // col 0 = 0E, increasing eastward

// Hourly snapshots available per day, in UTC.
inline constexpr std::array<int, 8> HOURS = {0, 3, 6, 9, 12, 15, 18, 21};

// Minimum number of valid hourly timesteps (per day, and per cell) required
// to produce a non-NaN output. Overridable via --min-valid-hours.
inline constexpr int DEFAULT_MIN_VALID_HOURS = 2;

} // namespace weather::config
