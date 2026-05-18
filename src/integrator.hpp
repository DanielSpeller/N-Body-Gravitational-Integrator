#pragma once

#include "body.hpp"
#include "physics.hpp"
#include <vector>
#include <array>

// Advance the system by one time step dt using velocity-Verlet integration.
//
// This is a symplectic second-order method that conserves a shadow Hamiltonian,
// giving bounded (oscillating) energy error rather than the secular drift of
// non-symplectic methods like RK4.
//
// Algorithm:
//   a(t)         = acc_current  (supplied by caller)
//   v(t+dt/2)    = v(t) + a(t) * dt/2
//   r(t+dt)      = r(t) + v(t+dt/2) * dt
//   a(t+dt)      = compute_forces(r(t+dt))
//   v(t+dt)      = v(t+dt/2) + a(t+dt) * dt/2
//
// acc_current is updated in-place to hold a(t+dt), so the next call reuses it
// without an extra force evaluation.
inline void verlet_step(
    std::vector<Body>&                  bodies,
    std::vector<std::array<double, 3>>& acc_current,
    double                              dt,
    double                              G)
{
    const std::size_t n = bodies.size();

    // Half-kick: v(t) -> v(t + dt/2)
    for (std::size_t i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k)
            bodies[i].vel[k] += acc_current[i][k] * 0.5 * dt;

    // Full drift: r(t) -> r(t + dt)
    for (std::size_t i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k)
            bodies[i].pos[k] += bodies[i].vel[k] * dt;

    // Recompute forces at new positions
    acc_current = compute_forces(bodies, G);

    // Half-kick: v(t + dt/2) -> v(t + dt)
    for (std::size_t i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k)
            bodies[i].vel[k] += acc_current[i][k] * 0.5 * dt;
}
