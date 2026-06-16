#include "io.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

#include <cnpy.h>

namespace weather::io {

namespace fs = std::filesystem;

namespace {

using RowMajorMatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

std::string zero_pad(int value, int width) {
    std::ostringstream oss;
    oss << std::setw(width) << std::setfill('0') << value;
    return oss.str();
}

// cnpy::parse_npy_header discards the numpy dtype "kind" character (it only
// keeps the byte width), so a 1-byte field cannot be told apart as int8 vs
// uint8, nor a 2-byte field as float16 vs int16. The real GFS-Wave archive
// mixes float16 continuous fields with int8-coded compass directions in the
// same hour directory, so the kind character has to be recovered directly
// from the on-disk header dict (e.g. "'descr': '<f2'" or "'descr': '|i1'").
struct DtypeInfo {
    char kind;       // 'f' (float), 'i' (signed int), 'u' (unsigned int)
    size_t word_size; // bytes per element
    std::string descr; // raw numpy dtype string, for error messages
};

DtypeInfo read_dtype_info(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        throw std::runtime_error("weather::io: cannot open " + path);
    }

    char magic[6];
    f.read(magic, sizeof(magic));
    char major = 0, minor = 0;
    f.read(&major, 1);
    f.read(&minor, 1);
    uint16_t header_len = 0;
    f.read(reinterpret_cast<char*>(&header_len), sizeof(header_len));
    std::string header(header_len, '\0');
    f.read(header.data(), header_len);
    if (!f) {
        throw std::runtime_error("weather::io: failed to read .npy header in " + path);
    }

    const size_t key = header.find("descr");
    if (key == std::string::npos) {
        throw std::runtime_error("weather::io: missing 'descr' field in " + path);
    }
    // header looks like: ... 'descr': '<f2', ... -- skip the key's own
    // closing quote, then take the value between the next pair of quotes.
    const size_t key_close_quote = header.find('\'', key + 5);
    const size_t value_open_quote = header.find('\'', key_close_quote + 1);
    const size_t value_close_quote = header.find('\'', value_open_quote + 1);
    if (value_open_quote == std::string::npos || value_close_quote == std::string::npos) {
        throw std::runtime_error("weather::io: malformed 'descr' field in " + path);
    }
    const std::string descr =
        header.substr(value_open_quote + 1, value_close_quote - value_open_quote - 1);

    if (descr.size() < 2) {
        throw std::runtime_error("weather::io: malformed dtype string '" + descr +
                                  "' in " + path);
    }
    // byte-order prefix is one of '<', '>', '|', '='
    const char kind = descr[1];
    const size_t word_size = static_cast<size_t>(std::stoi(descr.substr(2)));
    return {kind, word_size, descr};
}

// IEEE 754 binary16 -> double, since neither AppleClang's nor every target
// GCC's <stdfloat> reliably exposes std::float16_t.
double half_to_double(uint16_t bits) {
    const uint16_t sign = (bits >> 15) & 0x1;
    const uint16_t exponent = (bits >> 10) & 0x1F;
    const uint16_t mantissa = bits & 0x3FF;

    double value;
    if (exponent == 0) {
        value = std::ldexp(static_cast<double>(mantissa), -24); // subnormal
    } else if (exponent == 0x1F) {
        value = (mantissa == 0) ? std::numeric_limits<double>::infinity()
                                 : std::numeric_limits<double>::quiet_NaN();
    } else {
        value = std::ldexp(1.0 + static_cast<double>(mantissa) / 1024.0, exponent - 15);
    }
    return sign ? -value : value;
}

template <typename T>
Eigen::MatrixXd cast_buffer_to_double(const cnpy::NpyArray& arr) {
    const T* raw = arr.data<T>();
    const Eigen::Index rows = static_cast<Eigen::Index>(arr.shape[0]);
    const Eigen::Index cols = static_cast<Eigen::Index>(arr.shape[1]);
    RowMajorMatrixXd row_major(rows, cols);
    for (Eigen::Index idx = 0; idx < rows * cols; ++idx) {
        row_major.data()[idx] = static_cast<double>(raw[idx]);
    }
    return Eigen::MatrixXd(row_major);
}

Eigen::MatrixXd cast_half_buffer_to_double(const cnpy::NpyArray& arr) {
    const uint16_t* raw = arr.data<uint16_t>();
    const Eigen::Index rows = static_cast<Eigen::Index>(arr.shape[0]);
    const Eigen::Index cols = static_cast<Eigen::Index>(arr.shape[1]);
    RowMajorMatrixXd row_major(rows, cols);
    for (Eigen::Index idx = 0; idx < rows * cols; ++idx) {
        row_major.data()[idx] = half_to_double(raw[idx]);
    }
    return Eigen::MatrixXd(row_major);
}

} // namespace

Eigen::MatrixXd load_npy(const std::string& path) {
    if (!fs::exists(path)) {
        throw std::runtime_error("weather::io::load_npy: file does not exist: " + path);
    }

    const DtypeInfo dtype = read_dtype_info(path);
    cnpy::NpyArray arr = cnpy::npy_load(path);

    if (arr.shape.size() != 2 || arr.shape[0] == 0 || arr.shape[1] == 0) {
        throw std::runtime_error(
            "weather::io::load_npy: expected a non-empty 2D grid in " + path);
    }

    // The GFS-Wave archive stores continuous fields as float16/float64 and
    // direction fields as int8/uint8 compass codes; upcast everything to
    // double immediately so all downstream math is double-precision, per
    // AGENTS.md.
    if (dtype.kind == 'f' && dtype.word_size == 8) return cast_buffer_to_double<double>(arr);
    if (dtype.kind == 'f' && dtype.word_size == 4) return cast_buffer_to_double<float>(arr);
    if (dtype.kind == 'f' && dtype.word_size == 2) return cast_half_buffer_to_double(arr);
    if (dtype.kind == 'i' && dtype.word_size == 1) return cast_buffer_to_double<int8_t>(arr);
    if (dtype.kind == 'u' && dtype.word_size == 1) return cast_buffer_to_double<uint8_t>(arr);

    throw std::runtime_error("weather::io::load_npy: unsupported dtype '" + dtype.descr +
                              "' in " + path);
}

void save_npy(const std::string& path, const Eigen::MatrixXd& field) {
    fs::create_directories(fs::path(path).parent_path());

    // numpy's .npy format stores data in C (row-major) order by default;
    // convert before handing the raw buffer to cnpy::npy_save.
    RowMajorMatrixXd row_major = field;

    std::vector<size_t> shape = {static_cast<size_t>(field.rows()),
                                  static_cast<size_t>(field.cols())};
    cnpy::npy_save(path, row_major.data(), shape);
}

std::string path_for(const std::string& root, int year, int month, int day,
                      int hour, const std::string& varname) {
    fs::path p = fs::path(root) / zero_pad(year, 4) / zero_pad(month, 2) /
                 zero_pad(day, 2) / zero_pad(hour, 2) / (varname + ".npy");
    return p.string();
}

std::string path_for(const std::string& root, int year, int month, int day) {
    fs::path p =
        fs::path(root) / zero_pad(year, 4) / zero_pad(month, 2) / zero_pad(day, 2);
    return p.string();
}

} // namespace weather::io
