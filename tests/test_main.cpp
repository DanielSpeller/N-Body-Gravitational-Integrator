#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "body.hpp"
#include "physics.hpp"
#include "integrator.hpp"

#include <cmath>
#include <vector>

static constexpr double PI   = std::acos(-1.0);
static constexpr double G_AU = 4.0 * PI * PI;   // AU^3 / yr^2 / M_sun

// ── Helpers ──────────────────────────────────────────────────────────────────

static std::vector<Body> two_bodies_at_rest()
{
    Body b1{"A", 1.0, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    Body b2{"B", 1.0, {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    return {b1, b2};
}

// ═══════════════════════════════════════════════════════════════════════════
// Pairwise force / acceleration
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("compute_forces: two equal masses 1 AU apart")
{
    auto bodies = two_bodies_at_rest();
    auto acc    = compute_forces(bodies, G_AU);

    REQUIRE(acc.size() == 2);

    // Body 0 accelerates toward body 1 (+X); magnitude = G*m/r^2 = G_AU
    CHECK(acc[0][0] == doctest::Approx( G_AU).epsilon(1e-12));
    CHECK(acc[0][1] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(acc[0][2] == doctest::Approx(0.0).epsilon(1e-12));

    // Body 1 accelerates toward body 0 (-X) — Newton's 3rd law
    CHECK(acc[1][0] == doctest::Approx(-G_AU).epsilon(1e-12));
    CHECK(acc[1][1] == doctest::Approx(0.0).epsilon(1e-12));
    CHECK(acc[1][2] == doctest::Approx(0.0).epsilon(1e-12));
}

TEST_CASE("compute_forces: total momentum is conserved (sum m*a == 0)")
{
    // Equilateral triangle, equal masses — by symmetry net force = 0 system-wide
    Body b0{"A", 1.0, {0.0,            0.0,           0.0}, {0.0,0.0,0.0}};
    Body b1{"B", 1.0, {1.0,            0.0,           0.0}, {0.0,0.0,0.0}};
    Body b2{"C", 1.0, {0.5, std::sqrt(3.0)/2.0,       0.0}, {0.0,0.0,0.0}};
    std::vector<Body> bodies = {b0, b1, b2};

    auto acc = compute_forces(bodies, G_AU);

    // Z component must be zero (2D geometry)
    for (const auto& a : acc)
        CHECK(a[2] == doctest::Approx(0.0).epsilon(1e-12));

    // Sum of (m * a) must be zero in each direction
    double sum_x = 0.0, sum_y = 0.0;
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        sum_x += bodies[i].mass * acc[i][0];
        sum_y += bodies[i].mass * acc[i][1];
    }
    CHECK(sum_x == doctest::Approx(0.0).epsilon(1e-10));
    CHECK(sum_y == doctest::Approx(0.0).epsilon(1e-10));
}

// ═══════════════════════════════════════════════════════════════════════════
// Energy computation
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("total_energy: two bodies at rest — pure potential energy")
{
    auto bodies = two_bodies_at_rest();
    double E    = total_energy(bodies, G_AU);

    // KE = 0, PE = -G * m1 * m2 / r = -G_AU
    CHECK(E == doctest::Approx(-G_AU).epsilon(1e-12));
}

TEST_CASE("total_energy: single moving body — pure kinetic energy")
{
    Body b{"A", 2.0, {0.0, 0.0, 0.0}, {3.0, 4.0, 0.0}};
    std::vector<Body> bodies = {b};

    double E = total_energy(bodies, G_AU);
    // KE = 0.5 * 2 * (9 + 16) = 25
    CHECK(E == doctest::Approx(25.0).epsilon(1e-12));
}

// ═══════════════════════════════════════════════════════════════════════════
// Velocity-Verlet integrator
// ═══════════════════════════════════════════════════════════════════════════

TEST_CASE("verlet_step: free particle drifts by v*dt exactly")
{
    // Single body, no other bodies → zero acceleration
    Body b{"A", 1.0, {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
    std::vector<Body> bodies = {b};

    auto acc = compute_forces(bodies, G_AU);
    verlet_step(bodies, acc, 0.1, G_AU);

    CHECK(bodies[0].pos[0] == doctest::Approx(0.1).epsilon(1e-14));
    CHECK(bodies[0].vel[0] == doctest::Approx(1.0).epsilon(1e-14));
}

TEST_CASE("verlet_step: energy bounded over 1000 circular-orbit steps")
{
    // Earth-like orbit at 1 AU, v = 2*pi AU/yr
    Body sun  {"Sun",   1.0,   {0.0, 0.0, 0.0}, {0.0,       0.0, 0.0}};
    Body earth{"Earth", 3.0e-6,{1.0, 0.0, 0.0}, {0.0, 2.0*PI, 0.0}};
    std::vector<Body> bodies = {sun, earth};

    double E0  = total_energy(bodies, G_AU);
    auto   acc = compute_forces(bodies, G_AU);

    for (int i = 0; i < 1000; ++i)
        verlet_step(bodies, acc, 0.001, G_AU);

    double E1     = total_energy(bodies, G_AU);
    double dE_rel = std::abs((E1 - E0) / E0);

    // 1 year at dt=0.001 yr: relative energy error should be well below 1e-5
    CHECK(dE_rel < 1.0e-5);
}

TEST_CASE("verlet_step: angular momentum conserved over 1000 steps")
{
    Body sun  {"Sun",   1.0,   {0.0, 0.0, 0.0}, {0.0,       0.0, 0.0}};
    Body earth{"Earth", 3.0e-6,{1.0, 0.0, 0.0}, {0.0, 2.0*PI, 0.0}};
    std::vector<Body> bodies = {sun, earth};

    auto L0  = total_angular_momentum(bodies);
    auto acc = compute_forces(bodies, G_AU);

    for (int i = 0; i < 1000; ++i)
        verlet_step(bodies, acc, 0.001, G_AU);

    auto L1   = total_angular_momentum(bodies);
    double dLz = std::abs(L1[2] - L0[2]);

    CHECK(dLz < 1.0e-13);
}
