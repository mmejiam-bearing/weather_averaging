#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

#include "averaging.hpp"
#include "config.hpp"
#include "io.hpp"
#include "types.hpp"

namespace fs = std::filesystem;
using namespace weather;

namespace {

struct CliArgs {
    std::string data_dir;
    std::string output_dir;
    std::vector<int> years;
    std::vector<int> months;
    int threads = 0; // 0 means "leave OMP_NUM_THREADS as-is"
    int min_valid_hours = config::DEFAULT_MIN_VALID_HOURS;
};

[[noreturn]] void usage_error(const std::string& message) {
    std::cerr << "error: " << message << "\n\n"
              << "usage: weather_daily_avg \\\n"
              << "  --data-dir   <path> \\\n"
              << "  --output-dir <path> \\\n"
              << "  --years      <y1> [y2 ...] \\\n"
              << "  --months     <m1> [m2 ...] \\\n"
              << "  [--threads   <n>] \\\n"
              << "  [--min-valid-hours <n>]\n";
    std::exit(EXIT_FAILURE);
}

CliArgs parse_args(int argc, char** argv) {
    CliArgs args;
    std::vector<std::string> tokens(argv + 1, argv + argc);

    size_t i = 0;
    auto next_is_flag = [&](size_t idx) {
        return idx >= tokens.size() || (tokens[idx].size() > 1 && tokens[idx][0] == '-' &&
                                         !std::isdigit(static_cast<unsigned char>(tokens[idx][1])));
    };

    while (i < tokens.size()) {
        const std::string& flag = tokens[i];
        if (flag == "--data-dir") {
            if (++i >= tokens.size()) usage_error("--data-dir requires a value");
            args.data_dir = tokens[i++];
        } else if (flag == "--output-dir") {
            if (++i >= tokens.size()) usage_error("--output-dir requires a value");
            args.output_dir = tokens[i++];
        } else if (flag == "--years") {
            ++i;
            if (next_is_flag(i)) usage_error("--years requires at least one value");
            while (!next_is_flag(i)) args.years.push_back(std::stoi(tokens[i++]));
        } else if (flag == "--months") {
            ++i;
            if (next_is_flag(i)) usage_error("--months requires at least one value");
            while (!next_is_flag(i)) args.months.push_back(std::stoi(tokens[i++]));
        } else if (flag == "--threads") {
            if (++i >= tokens.size()) usage_error("--threads requires a value");
            args.threads = std::stoi(tokens[i++]);
        } else if (flag == "--min-valid-hours") {
            if (++i >= tokens.size()) usage_error("--min-valid-hours requires a value");
            args.min_valid_hours = std::stoi(tokens[i++]);
        } else {
            usage_error("unrecognized argument: " + flag);
        }
    }

    if (args.data_dir.empty()) usage_error("--data-dir is required");
    if (args.output_dir.empty()) usage_error("--output-dir is required");
    if (args.years.empty()) usage_error("--years requires at least one value");
    if (args.months.empty()) usage_error("--months requires at least one value");

    return args;
}

bool is_leap_year(int year) {
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

int days_in_month(int year, int month) {
    static constexpr int kDays[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (month < 1 || month > 12) return 0;
    int days = kDays[month - 1];
    if (month == 2 && is_leap_year(year)) days = 29;
    return days;
}

// Loads the WindSeaResult for a single hour, returning std::nullopt if any
// of the required input files are missing for that hour.
std::optional<WindSeaResult> load_hour(const std::string& data_dir, int year,
                                        int month, int day, int hour) {
    static constexpr const char* kVars[8] = {
        "sigwh", "wsh", "wsp", "wsd", "was", "wad", "pwh", "pwd"};
    for (const char* var : kVars) {
        if (!fs::exists(io::path_for(data_dir, year, month, day, hour, var))) {
            return std::nullopt;
        }
    }

    // wsd/wad/pwd are stored as int8 16-point compass codes in the archive
    // (wsd/pwd: 0-15 North-origin, wad: 1-16 North-origin), not raw degrees.
    WindSeaResult result;
    result.sigwh = io::load_npy(io::path_for(data_dir, year, month, day, hour, "sigwh"));
    result.wsh = io::load_npy(io::path_for(data_dir, year, month, day, hour, "wsh"));
    result.wsp = io::load_npy(io::path_for(data_dir, year, month, day, hour, "wsp"));
    result.wsd = compass_code_to_degrees(
        io::load_npy(io::path_for(data_dir, year, month, day, hour, "wsd")), 0.0);
    result.was = io::load_npy(io::path_for(data_dir, year, month, day, hour, "was"));
    result.wad = compass_code_to_degrees(
        io::load_npy(io::path_for(data_dir, year, month, day, hour, "wad")), 1.0);
    result.pwh = io::load_npy(io::path_for(data_dir, year, month, day, hour, "pwh"));
    result.pwd = compass_code_to_degrees(
        io::load_npy(io::path_for(data_dir, year, month, day, hour, "pwd")), 0.0);
    return result;
}

void process_day(const CliArgs& args, int year, int month, int day) {
    std::vector<WindSeaResult> hourly_snapshots;
    for (int hour : config::HOURS) {
        std::optional<WindSeaResult> snap = load_hour(args.data_dir, year, month, day, hour);
        if (snap) hourly_snapshots.push_back(std::move(*snap));
    }

    if (static_cast<int>(hourly_snapshots.size()) < args.min_valid_hours) {
        if (!hourly_snapshots.empty()) {
            #pragma omp critical(stderr_io)
            {
                std::cerr << "warning: " << year << "-" << month << "-" << day
                          << ": only " << hourly_snapshots.size()
                          << " valid hourly timestep(s) (< " << args.min_valid_hours
                          << " required); skipping day\n";
            }
        }
        return;
    }

    DailyAverages avg = daily_energy_average(hourly_snapshots, args.min_valid_hours);

    const std::string out_dir = io::path_for(args.output_dir, year, month, day);
    io::save_npy(out_dir + "/sigwh.npy", avg.sigwh);
    io::save_npy(out_dir + "/wsh.npy", avg.wsh);
    io::save_npy(out_dir + "/wsp.npy", avg.wsp);
    io::save_npy(out_dir + "/wsd.npy", avg.wsd);
    io::save_npy(out_dir + "/was.npy", avg.was);
    io::save_npy(out_dir + "/wad.npy", avg.wad);
    io::save_npy(out_dir + "/pwd.npy", avg.pwd);
    io::save_npy(out_dir + "/swell_residual.npy", avg.swell_residual);
}

} // namespace

int main(int argc, char** argv) {
    CliArgs args = parse_args(argc, argv);

    if (args.threads > 0) {
        omp_set_num_threads(args.threads);
    }

    struct DayTask {
        int year;
        int month;
        int day;
    };

    std::vector<DayTask> tasks;
    for (int year : args.years) {
        for (int month : args.months) {
            const int nd = days_in_month(year, month);
            for (int day = 1; day <= nd; ++day) tasks.push_back({year, month, day});
        }
    }

    #pragma omp parallel for schedule(dynamic)
    for (size_t idx = 0; idx < tasks.size(); ++idx) {
        const DayTask& t = tasks[idx];
        try {
            process_day(args, t.year, t.month, t.day);
        } catch (const std::exception& e) {
            #pragma omp critical(stderr_io)
            {
                std::cerr << "error: " << t.year << "-" << t.month << "-" << t.day
                          << ": " << e.what() << "\n";
            }
        }
    }

    return EXIT_SUCCESS;
}
