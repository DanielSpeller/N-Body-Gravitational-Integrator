#pragma once

#include "body.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <array>

// Read bodies from a whitespace-delimited file.
// Format per non-comment line:
//   name  mass  x  y  z  vx  vy  vz
// Lines beginning with '#' (after optional whitespace) are ignored.
inline std::vector<Body> read_bodies_from_file(const std::string& path)
{
    std::ifstream file(path);
    if (!file.is_open())
        throw std::runtime_error("Cannot open initial conditions file: " + path);

    std::vector<Body> bodies;
    std::string line;

    while (std::getline(file, line)) {
        std::size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos || line[first] == '#')
            continue;

        std::istringstream ss(line);
        Body b;
        if (!(ss >> b.name >> b.mass
                 >> b.pos[0] >> b.pos[1] >> b.pos[2]
                 >> b.vel[0] >> b.vel[1] >> b.vel[2])) {
            throw std::runtime_error("Malformed line in initial conditions: " + line);
        }
        bodies.push_back(b);
    }

    if (bodies.empty())
        throw std::runtime_error("No bodies found in file: " + path);

    return bodies;
}

// Write the CSV header row.
// Columns: time, name_x, name_y, name_z (per body), energy, L_x, L_y, L_z
inline void write_csv_header(std::ostream& out, const std::vector<Body>& bodies)
{
    out << "time";
    for (const auto& b : bodies)
        out << "," << b.name << "_x"
            << "," << b.name << "_y"
            << "," << b.name << "_z";
    out << ",energy,L_x,L_y,L_z\n";
}

// Write one CSV data row at the given simulation time.
inline void write_csv_row(std::ostream&               out,
                          double                      t,
                          const std::vector<Body>&    bodies,
                          double                      energy,
                          const std::array<double, 3>& L)
{
    out << t;
    for (const auto& b : bodies)
        out << "," << b.pos[0] << "," << b.pos[1] << "," << b.pos[2];
    out << "," << energy
        << "," << L[0] << "," << L[1] << "," << L[2]
        << "\n";
}
