#include "body.hpp"
#include "physics.hpp"
#include "integrator.hpp"
#include "io.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <filesystem>
#include <stdexcept>

static constexpr double PI         = std::acos(-1.0);
static constexpr double G_DEFAULT  = 4.0 * PI * PI;   // AU^3 / yr^2 / M_sun

struct Config {
    std::string input_file;
    double      dt     = 0.001;
    double      tmax   = 100.0;
    double      G      = G_DEFAULT;
    std::string output = "output/trajectory.csv";
    int         stride = 10;
};

static Config parse_args(int argc, char* argv[])
{
    Config cfg;
    if (argc < 2) {
        std::cerr << "Usage: nbody <initial_conditions.txt> "
                     "[--dt yr] [--tmax yr] [--G val] [--out file.csv] [--stride N]\n";
        std::exit(1);
    }
    cfg.input_file = argv[1];
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if      (arg == "--dt"     && i+1 < argc) { cfg.dt     = std::stod(argv[++i]); }
        else if (arg == "--tmax"   && i+1 < argc) { cfg.tmax   = std::stod(argv[++i]); }
        else if (arg == "--G"      && i+1 < argc) { cfg.G      = std::stod(argv[++i]); }
        else if (arg == "--out"    && i+1 < argc) { cfg.output = argv[++i]; }
        else if (arg == "--stride" && i+1 < argc) { cfg.stride = std::stoi(argv[++i]); }
        else {
            std::cerr << "Unknown argument: " << arg << "\n";
            std::exit(1);
        }
    }
    return cfg;
}

int main(int argc, char* argv[])
{
    const Config cfg = parse_args(argc, argv);

    std::vector<Body> bodies;
    try {
        bodies = read_bodies_from_file(cfg.input_file);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Loaded " << bodies.size() << " bodies from " << cfg.input_file << "\n";
    std::cout << "G = " << cfg.G << "  dt = " << cfg.dt
              << " yr  tmax = " << cfg.tmax << " yr\n";

    // Create output directory if needed
    {
        std::filesystem::path out_path(cfg.output);
        if (out_path.has_parent_path())
            std::filesystem::create_directories(out_path.parent_path());
    }

    std::ofstream csv(cfg.output);
    if (!csv.is_open()) {
        std::cerr << "Error: cannot open output file: " << cfg.output << "\n";
        return 1;
    }
    csv.precision(15);
    write_csv_header(csv, bodies);

    // Initial state
    auto   acc = compute_forces(bodies, cfg.G);
    double E0  = total_energy(bodies, cfg.G);
    auto   L0  = total_angular_momentum(bodies);

    std::cout << "Initial energy:   " << E0 << "\n";
    std::cout << "Initial |L|:      "
              << std::sqrt(L0[0]*L0[0] + L0[1]*L0[1] + L0[2]*L0[2]) << "\n";

    write_csv_row(csv, 0.0, bodies, E0, L0);

    const long long total_steps    = static_cast<long long>(cfg.tmax / cfg.dt);
    const long long report_every   = (total_steps / 10 < 1) ? 1 : total_steps / 10;

    for (long long step = 1; step <= total_steps; ++step) {
        verlet_step(bodies, acc, cfg.dt, cfg.G);

        if (step % cfg.stride == 0) {
            double t = step * cfg.dt;
            double E = total_energy(bodies, cfg.G);
            auto   L = total_angular_momentum(bodies);
            write_csv_row(csv, t, bodies, E, L);

            if (step % (report_every * cfg.stride) == 0) {
                double dE = (E - E0) / std::abs(E0);
                std::cout << "t=" << t << " yr  dE/E0=" << dE << "\n";
            }
        }
    }

    double E_final = total_energy(bodies, cfg.G);
    std::cout << "Final energy:     " << E_final << "\n";
    std::cout << "Relative drift:   " << (E_final - E0) / std::abs(E0) << "\n";
    std::cout << "Output written to: " << cfg.output << "\n";

    return 0;
}
