#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "averaging.hpp"
#include "types.hpp"

using namespace weather;

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kEps = 1e-9;

// Builds a single-cell (1x1) EnergyStack from a list of hourly values.
EnergyStack make_stack(std::initializer_list<double> values) {
    EnergyStack stack;
    for (double v : values) {
        Eigen::MatrixXd cell(1, 1);
        cell(0, 0) = v;
        stack.snapshots.push_back(cell);
    }
    return stack;
}

} // namespace

// --- grid geo-referencing constants (config.hpp) ---------------------------

TEST(GridGeoConstants, WaveGridLastRowReachesEmpiricallyDerivedNorthBound) {
    // Guards against a typo silently desynchronizing the empirically-derived
    // row-0 latitude from the row count -- see README.md Section 2. The
    // wave grid spans -75 (row 0) to +80 (last row).
    const double last_row_lat = config::WAVE_GRID_LAT_AT_ROW0 +
                                 (config::WAVE_GRID_ROWS - 1) * config::GRID_LAT_STEP_DEG;
    EXPECT_NEAR(last_row_lat, 80.0, kEps);
}

TEST(GridGeoConstants, WindGridSpansFullPoleToPoleRange) {
    // The wind grid has exactly enough rows for a full -90..90 span at 0.25
    // deg resolution with no truncation, unlike the wave grid.
    const double last_row_lat = config::WIND_GRID_LAT_AT_ROW0 +
                                 (config::WIND_GRID_ROWS - 1) * config::GRID_LAT_STEP_DEG;
    EXPECT_NEAR(config::WIND_GRID_LAT_AT_ROW0, -90.0, kEps);
    EXPECT_NEAR(last_row_lat, 90.0, kEps);
}

// --- energy_average -------------------------------------------------------

TEST(EnergyAverage, KnownValuesMatchAnalyticRms) {
    // README 4.1: H = [0, 4] -> H_eff = sqrt((0^2 + 4^2) / 2) = 2.828427...
    EnergyStack stack = make_stack({0.0, 4.0});
    Eigen::MatrixXd result = energy_average(stack, /*min_valid_hours=*/2);
    EXPECT_NEAR(result(0, 0), std::sqrt(8.0), kEps);
}

TEST(EnergyAverage, AllNaNInputProducesNaNOutput) {
    EnergyStack stack = make_stack({kNaN, kNaN, kNaN});
    Eigen::MatrixXd result = energy_average(stack, /*min_valid_hours=*/2);
    EXPECT_TRUE(std::isnan(result(0, 0)));
}

TEST(EnergyAverage, SingleValidTimestepIsValidWhenThresholdIsOne) {
    // Demonstrates the per-cell count>0 mechanism in isolation, independent
    // of the day-level min-valid-hours default used by the CLI.
    EnergyStack stack = make_stack({3.5, kNaN, kNaN});
    Eigen::MatrixXd result = energy_average(stack, /*min_valid_hours=*/1);
    EXPECT_NEAR(result(0, 0), 3.5, kEps);
}

TEST(EnergyAverage, BelowMinValidHoursThresholdProducesNaN) {
    // Section 5: fewer than min_valid_hours valid timesteps -> NaN,
    // regardless of the values present.
    EnergyStack stack = make_stack({3.5, kNaN, kNaN});
    Eigen::MatrixXd result = energy_average(stack, /*min_valid_hours=*/2);
    EXPECT_TRUE(std::isnan(result(0, 0)));
}

TEST(EnergyAverage, CellsAreIndependentAcrossAGrid) {
    EnergyStack stack;
    Eigen::MatrixXd h0(2, 2), h1(2, 2);
    h0 << 0.0, kNaN, 2.0, kNaN;
    h1 << 4.0, 5.0, 2.0, kNaN;
    stack.snapshots = {h0, h1};

    Eigen::MatrixXd result = energy_average(stack, /*min_valid_hours=*/1);
    EXPECT_NEAR(result(0, 0), std::sqrt((0.0 + 16.0) / 2.0), kEps);
    EXPECT_NEAR(result(0, 1), 5.0, kEps); // only one valid timestep
    EXPECT_NEAR(result(1, 0), std::sqrt((4.0 + 4.0) / 2.0), kEps);
    EXPECT_TRUE(std::isnan(result(1, 1))); // both timesteps NaN
}

// --- energy_weighted_period -------------------------------------------------

TEST(EnergyWeightedPeriod, ZeroHeightTimestepContributesZeroWeight) {
    // If a calm (H=0) timestep were given equal weight, the arithmetic mean
    // of [1000, 6] would dominate the result. The energy weighting must
    // exclude it entirely.
    EnergyStack heights = make_stack({0.0, 4.0});
    EnergyStack periods = make_stack({1000.0, 6.0});
    Eigen::MatrixXd result = energy_weighted_period(heights, periods, /*min_valid_hours=*/1);
    EXPECT_NEAR(result(0, 0), 6.0, kEps);
}

TEST(EnergyWeightedPeriod, AllCalmProducesNaN) {
    EnergyStack heights = make_stack({0.0, 0.0});
    EnergyStack periods = make_stack({8.0, 9.0});
    Eigen::MatrixXd result = energy_weighted_period(heights, periods, /*min_valid_hours=*/1);
    EXPECT_TRUE(std::isnan(result(0, 0)));
}

TEST(EnergyWeightedPeriod, PeriodNaNAtAValidHeightTimestepIsSkippedInNumerator) {
    EnergyStack heights = make_stack({2.0, 4.0});
    EnergyStack periods = make_stack({kNaN, 6.0});
    Eigen::MatrixXd result = energy_weighted_period(heights, periods, /*min_valid_hours=*/1);
    // weight sum = 4 + 16 = 20; numerator = 16*6 = 96 (first term skipped)
    EXPECT_NEAR(result(0, 0), 96.0 / 20.0, kEps);
}

// --- energy_weighted_circular_mean ------------------------------------------

TEST(EnergyWeightedCircularMean, NearNorthWrapAround) {
    // README 4.4: mean of 359 deg and 1 deg must be ~0 deg, not 180 deg.
    EnergyStack heights = make_stack({1.0, 1.0});
    EnergyStack angles = make_stack({359.0, 1.0});
    Eigen::MatrixXd result = energy_weighted_circular_mean(heights, angles, /*min_valid_hours=*/1);
    EXPECT_NEAR(result(0, 0), 0.0, 1e-6);
}

TEST(EnergyWeightedCircularMean, OppositeDirectionsCancelWithoutProducingNaN) {
    // README 4.4.1: when the mean resultant length R is ~0 (e.g. equally
    // weighted opposing directions), the documented behavior is to output
    // the raw atan2 result rather than NaN -- the low energy at such cells
    // is expected to be reflected by H_eff instead.
    EnergyStack heights = make_stack({1.0, 1.0});
    EnergyStack angles = make_stack({90.0, 270.0});
    Eigen::MatrixXd result = energy_weighted_circular_mean(heights, angles, /*min_valid_hours=*/1);
    EXPECT_FALSE(std::isnan(result(0, 0)));
}

TEST(EnergyWeightedCircularMean, HigherEnergyDirectionDominates) {
    EnergyStack heights = make_stack({1.0, 3.0});
    EnergyStack angles = make_stack({0.0, 90.0});
    Eigen::MatrixXd result = energy_weighted_circular_mean(heights, angles, /*min_valid_hours=*/1);
    // weight 9 at 90 deg vs weight 1 at 0 deg -> mean should lean strongly
    // toward 90 deg.
    EXPECT_GT(result(0, 0), 45.0);
}

TEST(EnergyWeightedCircularMean, JointNaNGuardSkipsTimestepEvenIfOtherFieldValid) {
    EnergyStack heights = make_stack({1.0, kNaN});
    EnergyStack angles = make_stack({kNaN, 45.0});
    Eigen::MatrixXd result = energy_weighted_circular_mean(heights, angles, /*min_valid_hours=*/1);
    EXPECT_TRUE(std::isnan(result(0, 0))); // no jointly-valid timestep at all
}

// --- swell_residual ----------------------------------------------------------

TEST(SwellResidual, PositiveResidualMatchesAnalyticFormula) {
    Eigen::MatrixXd sigwh(1, 1), wsh(1, 1);
    sigwh(0, 0) = 5.0;
    wsh(0, 0) = 3.0;
    Eigen::MatrixXd result = swell_residual(sigwh, wsh);
    EXPECT_NEAR(result(0, 0), std::sqrt(25.0 - 9.0), kEps);
}

TEST(SwellResidual, EpsilonNegativeRadicandClampsToZeroNotNaN) {
    Eigen::MatrixXd sigwh(1, 1), wsh(1, 1);
    sigwh(0, 0) = 3.0;
    wsh(0, 0) = 3.0 + 1e-10; // wsh > sigwh due to floating-point noise
    Eigen::MatrixXd result = swell_residual(sigwh, wsh);
    EXPECT_EQ(result(0, 0), 0.0);
}

TEST(SwellResidual, NaNInEitherInputPropagates) {
    Eigen::MatrixXd sigwh(1, 1), wsh(1, 1);
    sigwh(0, 0) = kNaN;
    wsh(0, 0) = 2.0;
    Eigen::MatrixXd result = swell_residual(sigwh, wsh);
    EXPECT_TRUE(std::isnan(result(0, 0)));
}

// --- compass_code_to_degrees ---------------------------------------------

TEST(CompassCodeToDegrees, ZeroOriginCodesDecodeAs16PointCompass) {
    // wsd archive convention: codes 0-15, North-origin (0 -> 0 deg).
    Eigen::MatrixXd codes(1, 4);
    codes << 0.0, 4.0, 8.0, 15.0;
    Eigen::MatrixXd degrees = compass_code_to_degrees(codes, /*origin_code=*/0.0);
    EXPECT_NEAR(degrees(0, 0), 0.0, kEps);
    EXPECT_NEAR(degrees(0, 1), 90.0, kEps);
    EXPECT_NEAR(degrees(0, 2), 180.0, kEps);
    EXPECT_NEAR(degrees(0, 3), 337.5, kEps);
}

TEST(CompassCodeToDegrees, OneOriginCodesDecodeAs16PointCompass) {
    // wad archive convention: codes 1-16, North-origin (1 -> 0 deg).
    Eigen::MatrixXd codes(1, 4);
    codes << 1.0, 5.0, 9.0, 16.0;
    Eigen::MatrixXd degrees = compass_code_to_degrees(codes, /*origin_code=*/1.0);
    EXPECT_NEAR(degrees(0, 0), 0.0, kEps);
    EXPECT_NEAR(degrees(0, 1), 90.0, kEps);
    EXPECT_NEAR(degrees(0, 2), 180.0, kEps);
    EXPECT_NEAR(degrees(0, 3), 337.5, kEps);
}

TEST(CompassCodeToDegrees, NaNPropagates) {
    Eigen::MatrixXd codes(1, 1);
    codes(0, 0) = kNaN;
    Eigen::MatrixXd degrees = compass_code_to_degrees(codes, /*origin_code=*/0.0);
    EXPECT_TRUE(std::isnan(degrees(0, 0)));
}

// --- cross-snapshot shape validation --------------------------------------

TEST(EnergyAverage, MismatchedSnapshotShapesThrow) {
    EnergyStack stack;
    stack.snapshots.push_back(Eigen::MatrixXd::Zero(2, 2));
    stack.snapshots.push_back(Eigen::MatrixXd::Zero(3, 3));
    EXPECT_THROW(energy_average(stack, /*min_valid_hours=*/1), std::invalid_argument);
}

// --- daily_energy_average (integration) --------------------------------------

TEST(DailyEnergyAverage, CombinesAllFieldsConsistently) {
    auto make_hour = [](double sigwh, double wsh, double wsp, double wsd, double was,
                         double wad, double pwh, double pwd) {
        WindSeaResult r;
        r.sigwh = Eigen::MatrixXd::Constant(1, 1, sigwh);
        r.wsh = Eigen::MatrixXd::Constant(1, 1, wsh);
        r.wsp = Eigen::MatrixXd::Constant(1, 1, wsp);
        r.wsd = Eigen::MatrixXd::Constant(1, 1, wsd);
        r.was = Eigen::MatrixXd::Constant(1, 1, was);
        r.wad = Eigen::MatrixXd::Constant(1, 1, wad);
        r.pwh = Eigen::MatrixXd::Constant(1, 1, pwh);
        r.pwd = Eigen::MatrixXd::Constant(1, 1, pwd);
        return r;
    };

    std::vector<WindSeaResult> hours = {
        make_hour(/*sigwh=*/0.0, /*wsh=*/0.0, /*wsp=*/10.0, /*wsd=*/10.0, /*was=*/0.0,
                  /*wad=*/99.0, /*pwh=*/0.0, /*pwd=*/99.0),
        make_hour(/*sigwh=*/4.0, /*wsh=*/3.0, /*wsp=*/8.0, /*wsd=*/20.0, /*was=*/6.0,
                  /*wad=*/45.0, /*pwh=*/5.0, /*pwd=*/70.0),
    };

    DailyAverages out = daily_energy_average(hours, /*min_valid_hours=*/2);

    EXPECT_NEAR(out.sigwh(0, 0), std::sqrt((0.0 + 16.0) / 2.0), kEps);
    EXPECT_NEAR(out.wsh(0, 0), std::sqrt((0.0 + 9.0) / 2.0), kEps);
    // wsp/wsd are weighted by wsh^2, not sigwh^2: only the second hour has
    // nonzero wsh, so it alone determines wsp and wsd.
    EXPECT_NEAR(out.wsp(0, 0), 8.0, kEps);
    EXPECT_NEAR(out.wsd(0, 0), 20.0, kEps);
    EXPECT_NEAR(out.was(0, 0), std::sqrt((0.0 + 36.0) / 2.0), kEps);
    // wad is weighted by was^2: only the second hour has nonzero was, so its
    // wad (45 deg) alone determines the result, ignoring the first hour's
    // wad value (99 deg) since its weight is zero.
    EXPECT_NEAR(out.wad(0, 0), 45.0, kEps);
    // pwd is weighted by pwh^2 (its own paired height), not wsh^2 or
    // sigwh^2: only the second hour has nonzero pwh, so its pwd (70 deg)
    // alone determines the result.
    EXPECT_NEAR(out.pwd(0, 0), 70.0, kEps);
    EXPECT_NEAR(out.swell_residual(0, 0),
                std::sqrt(out.sigwh(0, 0) * out.sigwh(0, 0) -
                          out.wsh(0, 0) * out.wsh(0, 0)),
                kEps);
}
