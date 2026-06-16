#pragma once

#include <vector>

#include <Eigen/Dense>

#include "config.hpp"

namespace weather {

// Dimensions of a regular lat/lon grid. The pipeline has no single grid
// shared by every field -- wave fields and wind fields are archived on
// different row counts (see config::WAVE_GRID_ROWS / WIND_GRID_ROWS).
struct GridShape {
    int rows;
    int cols;
};

inline constexpr GridShape WAVE_GRID{config::WAVE_GRID_ROWS, config::GRID_COLS};
inline constexpr GridShape WIND_GRID{config::WIND_GRID_ROWS, config::GRID_COLS};

// One hourly snapshot of all retained GFS-Wave input variables, loaded from
// a single <data-dir>/<Y>/<M>/<D>/<H>/ directory.
struct WindSeaResult {
    Eigen::MatrixXd sigwh; // significant wave height (combined sea+swell), m
    Eigen::MatrixXd wsh;   // wind sea significant wave height, m
    Eigen::MatrixXd wsp;   // wind sea mean wave period, s
    Eigen::MatrixXd wsd;   // wind sea mean wave direction, degrees
    Eigen::MatrixXd was;   // wind speed at sea surface, m/s
    Eigen::MatrixXd wad;   // wind direction at sea surface, degrees
    Eigen::MatrixXd pwh;   // peak wave height, m (loaded only to weight pwd; not itself output)
    Eigen::MatrixXd pwd;   // peak wave direction, degrees
};

// A day's worth of hourly snapshots for a single scalar field, in chronological
// order, used as the input to the energy-averaging functions.
struct EnergyStack {
    std::vector<Eigen::MatrixXd> snapshots;
};

// The 8 daily-averaged output fields written to <output-dir>/<Y>/<M>/<D>/.
struct DailyAverages {
    Eigen::MatrixXd sigwh;           // H_eff,sig
    Eigen::MatrixXd wsh;             // H_eff,ws
    Eigen::MatrixXd wsp;             // T_eff,ws
    Eigen::MatrixXd wsd;             // theta_eff,ws
    Eigen::MatrixXd was;             // U_eff
    Eigen::MatrixXd wad;             // theta_eff,wind -- energy-weighted by was^2
    Eigen::MatrixXd pwd;             // theta_eff,peak -- energy-weighted by pwh^2
    Eigen::MatrixXd swell_residual;  // H_swell
};

} // namespace weather
