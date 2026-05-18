#pragma once

#include "body.hpp"
#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>

// Returns acceleration vectors for each body (length == bodies.size()).
// Uses Newton's 3rd law optimisation: each pair computed once.
inline std::vector<std::array<double, 3>>
compute_forces(const std::vector<Body>& bodies, double G)
{
    const std::size_t n = bodies.size();
    std::vector<std::array<double, 3>> acc(n, {0.0, 0.0, 0.0});

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            double dx = bodies[j].pos[0] - bodies[i].pos[0];
            double dy = bodies[j].pos[1] - bodies[i].pos[1];
            double dz = bodies[j].pos[2] - bodies[i].pos[2];

            double dist2 = dx*dx + dy*dy + dz*dz;
            if (dist2 == 0.0)
                throw std::runtime_error("Two bodies occupy the same position.");

            double dist3 = dist2 * std::sqrt(dist2);

            // a_i += G * m_j / |r_ij|^3 * (r_j - r_i)
            double ci = G * bodies[j].mass / dist3;
            acc[i][0] += ci * dx;
            acc[i][1] += ci * dy;
            acc[i][2] += ci * dz;

            // Newton's 3rd law
            double cj = G * bodies[i].mass / dist3;
            acc[j][0] -= cj * dx;
            acc[j][1] -= cj * dy;
            acc[j][2] -= cj * dz;
        }
    }
    return acc;
}

// Total mechanical energy: KE + PE
// E = sum_i (0.5 * m_i * |v_i|^2) - sum_{i<j} G*m_i*m_j / |r_ij|
inline double total_energy(const std::vector<Body>& bodies, double G)
{
    const std::size_t n = bodies.size();
    double E = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        double v2 = 0.0;
        for (int k = 0; k < 3; ++k)
            v2 += bodies[i].vel[k] * bodies[i].vel[k];
        E += 0.5 * bodies[i].mass * v2;

        for (std::size_t j = i + 1; j < n; ++j) {
            double dx = bodies[j].pos[0] - bodies[i].pos[0];
            double dy = bodies[j].pos[1] - bodies[i].pos[1];
            double dz = bodies[j].pos[2] - bodies[i].pos[2];
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            E -= G * bodies[i].mass * bodies[j].mass / dist;
        }
    }
    return E;
}

// Total angular momentum vector L = sum_i m_i * (r_i x v_i)
inline std::array<double, 3>
total_angular_momentum(const std::vector<Body>& bodies)
{
    std::array<double, 3> L = {0.0, 0.0, 0.0};
    for (const auto& b : bodies) {
        L[0] += b.mass * (b.pos[1]*b.vel[2] - b.pos[2]*b.vel[1]);
        L[1] += b.mass * (b.pos[2]*b.vel[0] - b.pos[0]*b.vel[2]);
        L[2] += b.mass * (b.pos[0]*b.vel[1] - b.pos[1]*b.vel[0]);
    }
    return L;
}
