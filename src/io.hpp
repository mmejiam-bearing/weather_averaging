#pragma once

#include <string>

#include <Eigen/Dense>

namespace weather::io {

// Loads a 2D .npy field from disk into a (col-major) Eigen::MatrixXd,
// upcasting to double and converting from cnpy's row-major buffer. Accepts
// float64, float32, float16, int8, and uint8 element types, since the
// GFS-Wave archive stores continuous fields as float16 and direction fields
// as int8/uint8 compass codes. Each field's own on-disk shape is treated as
// authoritative -- no fixed grid size is enforced here (e.g. wave fields and
// wind fields use different row counts in the real archive).
// Throws std::runtime_error if the file is missing, is not 2D, or uses an
// unsupported dtype.
Eigen::MatrixXd load_npy(const std::string& path);

// Writes a 721x1440 Eigen::MatrixXd field to disk as a row-major float64
// .npy file. Creates parent directories if they do not already exist.
void save_npy(const std::string& path, const Eigen::MatrixXd& field);

// Builds the path to a single hourly input variable file:
//   <root>/<YYYY>/<MM>/<DD>/<HH>/<varname>.npy
std::string path_for(const std::string& root, int year, int month, int day,
                      int hour, const std::string& varname);

// Builds the path to a daily output directory:
//   <root>/<YYYY>/<MM>/<DD>
std::string path_for(const std::string& root, int year, int month, int day);

} // namespace weather::io
