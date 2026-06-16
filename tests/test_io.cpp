#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#include <gtest/gtest.h>

#include <cnpy.h>

#include "config.hpp"
#include "io.hpp"

using namespace weather;
namespace fs = std::filesystem;

namespace {

fs::path make_temp_dir() {
    fs::path dir = fs::temp_directory_path() /
                    fs::path("weather_io_test_" + std::to_string(::getpid()));
    fs::create_directories(dir);
    return dir;
}

// Writes a minimal .npy v1.0 file with an arbitrary dtype string, for dtypes
// cnpy::npy_save<T> cannot express natively (e.g. float16, int32 fixtures
// used only to exercise weather::io::load_npy's dtype dispatch).
void write_raw_npy(const fs::path& path, const std::string& descr,
                    const std::vector<size_t>& shape,
                    const std::vector<char>& raw_bytes) {
    std::ostringstream dict;
    dict << "{'descr': '" << descr << "', 'fortran_order': False, 'shape': (";
    for (size_t i = 0; i < shape.size(); ++i) {
        if (i > 0) dict << ", ";
        dict << shape[i];
    }
    if (shape.size() == 1) dict << ",";
    dict << "), }";
    std::string dict_str = dict.str();
    size_t remainder = 16 - (10 + dict_str.size()) % 16;
    dict_str.append(remainder, ' ');
    dict_str.back() = '\n';

    std::ofstream f(path, std::ios::binary);
    f.put(static_cast<char>(0x93));
    f.write("NUMPY", 5);
    f.put(1);
    f.put(0);
    uint16_t header_len = static_cast<uint16_t>(dict_str.size());
    f.write(reinterpret_cast<const char*>(&header_len), sizeof(header_len));
    f.write(dict_str.data(), static_cast<std::streamsize>(dict_str.size()));
    f.write(raw_bytes.data(), static_cast<std::streamsize>(raw_bytes.size()));
}

template <typename T>
std::vector<char> to_bytes(const std::vector<T>& values) {
    std::vector<char> bytes(values.size() * sizeof(T));
    std::memcpy(bytes.data(), values.data(), bytes.size());
    return bytes;
}

} // namespace

TEST(IoPathFor, HourlyInputPathIsZeroPaddedAndStructured) {
    std::string path = io::path_for("/data", 2024, 1, 5, 3, "sigwh");
    EXPECT_EQ(path, "/data/2024/01/05/03/sigwh.npy");
}

TEST(IoPathFor, DailyOutputDirPathIsZeroPadded) {
    std::string path = io::path_for("/out", 2024, 11, 30);
    EXPECT_EQ(path, "/out/2024/11/30");
}

TEST(IoLoadNpy, MissingFileThrows) {
    EXPECT_THROW(io::load_npy("/nonexistent/path/sigwh.npy"), std::runtime_error);
}

TEST(IoLoadNpy, ArbitrarySmallGridShapeLoadsFine) {
    // Each field's own on-disk shape is authoritative; no fixed 721x1440
    // grid is enforced (the real archive mixes wave-grid and wind-grid
    // shapes within the same hour directory).
    fs::path dir = make_temp_dir();
    fs::path file = dir / "small_shape.npy";

    std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
    cnpy::npy_save(file.string(), data.data(), std::vector<size_t>{2, 2});

    Eigen::MatrixXd loaded = io::load_npy(file.string());
    ASSERT_EQ(loaded.rows(), 2);
    ASSERT_EQ(loaded.cols(), 2);
    EXPECT_DOUBLE_EQ(loaded(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(loaded(1, 1), 4.0);

    fs::remove_all(dir);
}

TEST(IoLoadNpy, OneDimensionalShapeThrows) {
    fs::path dir = make_temp_dir();
    fs::path file = dir / "one_d.npy";

    std::vector<double> data = {1.0, 2.0, 3.0, 4.0};
    cnpy::npy_save(file.string(), data);

    EXPECT_THROW(io::load_npy(file.string()), std::runtime_error);

    fs::remove_all(dir);
}

TEST(IoLoadNpy, UnsupportedDtypeThrows) {
    fs::path dir = make_temp_dir();
    fs::path file = dir / "int32.npy";

    std::vector<int32_t> data = {1, 2, 3, 4};
    write_raw_npy(file, "<i4", {2, 2}, to_bytes(data));

    EXPECT_THROW(io::load_npy(file.string()), std::runtime_error);

    fs::remove_all(dir);
}

TEST(IoLoadNpy, Int8CompassCodesUpcastToDouble) {
    fs::path dir = make_temp_dir();
    fs::path file = dir / "wsd.npy";

    std::vector<int8_t> codes = {-128, -3, 0, 5, 127};
    write_raw_npy(file, "|i1", {5, 1}, to_bytes(codes));

    Eigen::MatrixXd loaded = io::load_npy(file.string());
    ASSERT_EQ(loaded.rows(), 5);
    EXPECT_DOUBLE_EQ(loaded(0, 0), -128.0);
    EXPECT_DOUBLE_EQ(loaded(1, 0), -3.0);
    EXPECT_DOUBLE_EQ(loaded(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(loaded(3, 0), 5.0);
    EXPECT_DOUBLE_EQ(loaded(4, 0), 127.0);

    fs::remove_all(dir);
}

TEST(IoLoadNpy, UInt8CodesUpcastToDouble) {
    fs::path dir = make_temp_dir();
    fs::path file = dir / "uint8.npy";

    std::vector<uint8_t> codes = {0, 128, 255};
    write_raw_npy(file, "|u1", {3, 1}, to_bytes(codes));

    Eigen::MatrixXd loaded = io::load_npy(file.string());
    EXPECT_DOUBLE_EQ(loaded(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(loaded(1, 0), 128.0);
    EXPECT_DOUBLE_EQ(loaded(2, 0), 255.0);

    fs::remove_all(dir);
}

TEST(IoLoadNpy, Float16DecodesKnownBitPatterns) {
    // GFS-Wave continuous fields (sigwh, wsh, wsp, was) are stored as IEEE
    // binary16. Exercise zero, +/-normal, +/-infinity, NaN, and the smallest
    // subnormal against weather::io's manual half->double decoder.
    fs::path dir = make_temp_dir();
    fs::path file = dir / "float16.npy";

    std::vector<uint16_t> bits = {
        0x0000, // +0.0
        0x3C00, // 1.0
        0xC000, // -2.0
        0x7C00, // +Inf
        0xFC00, // -Inf
        0x7E00, // NaN
        0x0001, // smallest subnormal: 2^-24
    };
    write_raw_npy(file, "<f2", {7, 1}, to_bytes(bits));

    Eigen::MatrixXd loaded = io::load_npy(file.string());
    EXPECT_DOUBLE_EQ(loaded(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(loaded(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(loaded(2, 0), -2.0);
    EXPECT_TRUE(std::isinf(loaded(3, 0)) && loaded(3, 0) > 0);
    EXPECT_TRUE(std::isinf(loaded(4, 0)) && loaded(4, 0) < 0);
    EXPECT_TRUE(std::isnan(loaded(5, 0)));
    EXPECT_DOUBLE_EQ(loaded(6, 0), std::ldexp(1.0, -24));

    fs::remove_all(dir);
}

TEST(IoRoundTrip, SaveThenLoadPreservesValuesAndStorageOrder) {
    fs::path dir = make_temp_dir();
    fs::path file = dir / "roundtrip.npy";

    Eigen::MatrixXd field(config::WIND_GRID_ROWS, config::GRID_COLS);
    for (Eigen::Index i = 0; i < field.rows(); ++i) {
        for (Eigen::Index j = 0; j < field.cols(); ++j) {
            field(i, j) = static_cast<double>(i) * 1000.0 + static_cast<double>(j);
        }
    }
    field(0, 0) = std::numeric_limits<double>::quiet_NaN();

    io::save_npy(file.string(), field);
    Eigen::MatrixXd loaded = io::load_npy(file.string());

    ASSERT_EQ(loaded.rows(), field.rows());
    ASSERT_EQ(loaded.cols(), field.cols());
    EXPECT_TRUE(std::isnan(loaded(0, 0)));
    EXPECT_DOUBLE_EQ(loaded(5, 7), field(5, 7));
    EXPECT_DOUBLE_EQ(loaded(field.rows() - 1, field.cols() - 1),
                      field(field.rows() - 1, field.cols() - 1));

    fs::remove_all(dir);
}

TEST(IoSaveNpy, CreatesParentDirectoriesIfMissing) {
    fs::path dir = make_temp_dir();
    fs::path nested = dir / "nested" / "subdir" / "field.npy";

    Eigen::MatrixXd field = Eigen::MatrixXd::Constant(config::WIND_GRID_ROWS, config::GRID_COLS, 1.0);
    io::save_npy(nested.string(), field);

    EXPECT_TRUE(fs::exists(nested));

    fs::remove_all(dir);
}
