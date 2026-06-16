#include "averaging.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace weather {

namespace {

constexpr double PI = 3.14159265358979323846;

void require_uniform_shape(const EnergyStack& stack, const char* fn) {
    const Eigen::MatrixXd& first = stack.snapshots.front();
    for (const Eigen::MatrixXd& snap : stack.snapshots) {
        if (snap.rows() != first.rows() || snap.cols() != first.cols()) {
            throw std::invalid_argument(
                std::string(fn) + ": hourly snapshots have inconsistent shapes");
        }
    }
}

void require_nonempty_stack(const EnergyStack& stack, const char* fn) {
    if (stack.snapshots.empty()) {
        throw std::invalid_argument(std::string(fn) + ": empty stack");
    }
    require_uniform_shape(stack, fn);
}

void require_matching_stacks(const EnergyStack& a, const EnergyStack& b,
                              const char* fn) {
    if (a.snapshots.empty() || a.snapshots.size() != b.snapshots.size()) {
        throw std::invalid_argument(std::string(fn) +
                                     ": mismatched or empty stacks");
    }
    require_uniform_shape(a, fn);
    require_uniform_shape(b, fn);
    if (a.snapshots.front().rows() != b.snapshots.front().rows() ||
        a.snapshots.front().cols() != b.snapshots.front().cols()) {
        throw std::invalid_argument(std::string(fn) +
                                     ": height and secondary field shapes differ");
    }
}

} // namespace

Eigen::MatrixXd energy_average(const EnergyStack& stack, int min_valid_hours) {
    require_nonempty_stack(stack, "energy_average");

    const Eigen::Index rows = stack.snapshots.front().rows();
    const Eigen::Index cols = stack.snapshots.front().cols();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    Eigen::MatrixXd result(rows, cols);

    for (Eigen::Index i = 0; i < rows; ++i) {
        for (Eigen::Index j = 0; j < cols; ++j) {
            double sum_sq = 0.0;
            int count = 0;
            for (const Eigen::MatrixXd& snap : stack.snapshots) {
                const double h = snap(i, j);
                if (!std::isnan(h)) {
                    sum_sq += h * h;
                    ++count;
                }
            }
            result(i, j) = (count >= min_valid_hours)
                               ? std::sqrt(sum_sq / static_cast<double>(count))
                               : nan;
        }
    }
    return result;
}

Eigen::MatrixXd energy_weighted_period(const EnergyStack& height_stack,
                                        const EnergyStack& period_stack,
                                        int min_valid_hours) {
    require_matching_stacks(height_stack, period_stack, "energy_weighted_period");

    const Eigen::Index rows = height_stack.snapshots.front().rows();
    const Eigen::Index cols = height_stack.snapshots.front().cols();
    const size_t n = height_stack.snapshots.size();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    Eigen::MatrixXd result(rows, cols);

    for (Eigen::Index i = 0; i < rows; ++i) {
        for (Eigen::Index j = 0; j < cols; ++j) {
            double sum_sq = 0.0;
            double sum_Tw = 0.0;
            int count = 0;
            for (size_t t = 0; t < n; ++t) {
                const double h = height_stack.snapshots[t](i, j);
                if (std::isnan(h)) continue;
                const double w = h * h;
                sum_sq += w;
                ++count;
                const double T = period_stack.snapshots[t](i, j);
                if (!std::isnan(T)) sum_Tw += w * T;
            }
            result(i, j) = (count >= min_valid_hours && sum_sq > 0.0)
                               ? sum_Tw / sum_sq
                               : nan;
        }
    }
    return result;
}

Eigen::MatrixXd energy_weighted_circular_mean(const EnergyStack& height_stack,
                                               const EnergyStack& angle_deg_stack,
                                               int min_valid_hours) {
    require_matching_stacks(height_stack, angle_deg_stack,
                             "energy_weighted_circular_mean");

    const Eigen::Index rows = height_stack.snapshots.front().rows();
    const Eigen::Index cols = height_stack.snapshots.front().cols();
    const size_t n = height_stack.snapshots.size();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    Eigen::MatrixXd result(rows, cols);

    for (Eigen::Index i = 0; i < rows; ++i) {
        for (Eigen::Index j = 0; j < cols; ++j) {
            double sum_u = 0.0;
            double sum_v = 0.0;
            int count = 0;
            for (size_t t = 0; t < n; ++t) {
                const double h = height_stack.snapshots[t](i, j);
                const double theta = angle_deg_stack.snapshots[t](i, j);
                if (std::isnan(h) || std::isnan(theta)) continue;
                const double w = h * h;
                const double rad = theta * PI / 180.0;
                sum_u += w * std::sin(rad);
                sum_v += w * std::cos(rad);
                ++count;
            }
            if (count < min_valid_hours) {
                result(i, j) = nan;
            } else {
                // atan2 argument order is (sin, cos) per the meteorological
                // clockwise-from-North convention — do not swap to (y, x).
                result(i, j) = std::fmod(
                    std::atan2(sum_u, sum_v) * 180.0 / PI + 360.0, 360.0);
            }
        }
    }
    return result;
}

Eigen::MatrixXd swell_residual(const Eigen::MatrixXd& sigwh_daily,
                                const Eigen::MatrixXd& wsh_daily) {
    if (sigwh_daily.rows() != wsh_daily.rows() ||
        sigwh_daily.cols() != wsh_daily.cols()) {
        throw std::invalid_argument("swell_residual: mismatched field shapes");
    }

    const double nan = std::numeric_limits<double>::quiet_NaN();
    Eigen::MatrixXd result(sigwh_daily.rows(), sigwh_daily.cols());

    for (Eigen::Index i = 0; i < sigwh_daily.rows(); ++i) {
        for (Eigen::Index j = 0; j < sigwh_daily.cols(); ++j) {
            const double sig = sigwh_daily(i, j);
            const double ws = wsh_daily(i, j);
            if (std::isnan(sig) || std::isnan(ws)) {
                result(i, j) = nan;
            } else {
                result(i, j) = std::sqrt(std::max(sig * sig - ws * ws, 0.0));
            }
        }
    }
    return result;
}

DailyAverages daily_energy_average(
    const std::vector<WindSeaResult>& hourly_snapshots, int min_valid_hours) {
    if (hourly_snapshots.empty()) {
        throw std::invalid_argument("daily_energy_average: no hourly snapshots");
    }

    EnergyStack sigwh_stack, wsh_stack, wsp_stack, wsd_stack, was_stack, wad_stack,
        pwh_stack, pwd_stack;
    const size_t n = hourly_snapshots.size();
    sigwh_stack.snapshots.reserve(n);
    wsh_stack.snapshots.reserve(n);
    wsp_stack.snapshots.reserve(n);
    wsd_stack.snapshots.reserve(n);
    was_stack.snapshots.reserve(n);
    wad_stack.snapshots.reserve(n);
    pwh_stack.snapshots.reserve(n);
    pwd_stack.snapshots.reserve(n);

    for (const WindSeaResult& snap : hourly_snapshots) {
        sigwh_stack.snapshots.push_back(snap.sigwh);
        wsh_stack.snapshots.push_back(snap.wsh);
        wsp_stack.snapshots.push_back(snap.wsp);
        wsd_stack.snapshots.push_back(snap.wsd);
        was_stack.snapshots.push_back(snap.was);
        wad_stack.snapshots.push_back(snap.wad);
        pwh_stack.snapshots.push_back(snap.pwh);
        pwd_stack.snapshots.push_back(snap.pwd);
    }

    DailyAverages out;
    out.sigwh = energy_average(sigwh_stack, min_valid_hours);
    out.wsh = energy_average(wsh_stack, min_valid_hours);
    // wsp/wsd describe the wind-sea component specifically, so they are
    // weighted by wsh^2 (wind-sea height), not the combined sigwh^2.
    out.wsp = energy_weighted_period(wsh_stack, wsp_stack, min_valid_hours);
    out.wsd = energy_weighted_circular_mean(wsh_stack, wsd_stack, min_valid_hours);
    out.was = energy_average(was_stack, min_valid_hours);
    // wad is weighted by was^2, mirroring wsd's wsh^2 weighting -- consistent
    // with the wind-drag relation F_wind ~ U^2 (README Section 4.5).
    out.wad = energy_weighted_circular_mean(was_stack, wad_stack, min_valid_hours);
    // pwd is weighted by pwh^2 (its own paired peak-wave height), not sigwh^2,
    // for the same reason wsd is weighted by wsh^2 rather than sigwh^2: the
    // direction belongs to a specific spectral component, not the combined
    // sea state. pwh itself is not an output field.
    out.pwd = energy_weighted_circular_mean(pwh_stack, pwd_stack, min_valid_hours);
    out.swell_residual = swell_residual(out.sigwh, out.wsh);
    return out;
}

Eigen::MatrixXd compass_code_to_degrees(const Eigen::MatrixXd& codes,
                                         double origin_code) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    Eigen::MatrixXd result(codes.rows(), codes.cols());

    for (Eigen::Index i = 0; i < codes.rows(); ++i) {
        for (Eigen::Index j = 0; j < codes.cols(); ++j) {
            const double code = codes(i, j);
            if (std::isnan(code)) {
                result(i, j) = nan;
            } else {
                result(i, j) = std::fmod((code - origin_code) * 22.5 + 360.0, 360.0);
            }
        }
    }
    return result;
}

} // namespace weather
