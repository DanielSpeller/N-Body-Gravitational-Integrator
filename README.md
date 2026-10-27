# N-Body Gravitational Integrator

Simulates the gravitational motion of N point masses using the velocity-Verlet (leapfrog) integrator. Comes with initial conditions for the Sun-Earth-Jupiter system and Burrau's Pythagorean three-body problem.


---

## how to build it

You need CMake >= 3.16 and a C++17 compiler. MinGW (g++) works fine on Windows.

```bash
# configure (fetches the doctest header automatically)
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release

# build
cmake --build build
```

This produces `build/nbody.exe` and `build/nbody_tests.exe`.

---

## how to run it

```bash
# Sun-Earth-Jupiter: 100 years, dt=0.001 yr
./build/nbody.exe data/sun_earth_jupiter.txt --tmax 100 --dt 0.001 --out output/sej.csv

# Burrau's Pythagorean problem (G=1, dimensionless units)
./build/nbody.exe data/pythagorean.txt --G 1.0 --tmax 20 --dt 0.001 --out output/pyth.csv
```

All flags are optional:

| flag | default | description |
|---|---|---|
| `--dt` | 0.001 | time step in years |
| `--tmax` | 100.0 | end time in years |
| `--G` | 4π² | gravitational constant |
| `--out` | output/trajectory.csv | where to write the CSV |
| `--stride` | 10 | write every N-th step |

The `output/` directory is created automatically.

---

## running the tests

```bash
cd build && ctest --output-on-failure
```

Or run directly: `./build/nbody_tests.exe`

---

## what the output looks like

Each row in the CSV is one timestep:

```
time, Sun_x, Sun_y, Sun_z, Earth_x, ..., energy, L_x, L_y, L_z
```

You can plot trajectories and energy drift with:

```bash
python analysis/plot.py output/sej.csv
python analysis/plot.py output/pyth.csv --save   # writes PNG files
```

---

## the integrator

Velocity-Verlet (also called leapfrog) runs like this each step:

```
a(t)        = sum over pairs:  G * m_j * (r_j - r_i) / |r_j - r_i|^3
v(t+dt/2)   = v(t) + a(t) * dt/2         <- half-kick
r(t+dt)     = r(t) + v(t+dt/2) * dt      <- full drift
a(t+dt)     = recompute forces            <- one force eval per step
v(t+dt)     = v(t+dt/2) + a(t+dt) * dt/2 <- half-kick
```


---

## the initial conditions

**Sun-Earth-Jupiter** (`data/sun_earth_jupiter.txt`) — heliocentric circular orbit approximation, units of AU / yr / M_sun:

| body | mass (M_sun) | x (AU) | vy (AU/yr) |
|---|---|---|---|
| Sun | 1.0 | 0 | 0 |
| Earth | 3e-6 | 1.0 | 2*pi |
| Jupiter | 1e-3 | 5.2 | 2*pi/sqrt(5.2) |

**Pythagorean problem** (`data/pythagorean.txt`) — three masses at rest at the vertices of a 3-4-5 right triangle (Burrau 1913). The system is chaotic: a close triple encounter near t~15 ejects one body.

| body | mass | x | y |
|---|---|---|---|
| Body1 | 3 | 1 | 3 |
| Body2 | 4 | -2 | -1 |
| Body3 | 5 | 1 | -1 |

---

## files in this repo

```
.
├── CMakeLists.txt
├── src/
│   ├── body.hpp          Body struct (mass, position, velocity)
│   ├── physics.hpp       force computation, total energy, angular momentum
│   ├── integrator.hpp    velocity-Verlet step
│   ├── io.hpp            file I/O and CSV output
│   └── main.cpp          CLI driver
├── tests/
│   └── test_main.cpp     doctest unit tests
├── data/
│   ├── sun_earth_jupiter.txt
│   └── pythagorean.txt
└── analysis/
    └── plot.py           trajectory and energy drift plots
```

---

## references

- Burrau, C. (1913). Numerische Berechnung eines Spezialfalles des Dreikorper-Problems. *Astronomische Nachrichten*, 195(6), 113-118.
- Szebehely, V. & Peters, C.F. (1967). Complete Solution of a General Problem of Three Bodies. *Astronomical Journal*, 72, 876.
- Leimkuhler, B. & Reich, S. (2004). *Simulating Hamiltonian Dynamics*. Cambridge.
- Hairer, E., Lubich, C. & Wanner, G. (2006). *Geometric Numerical Integration*. Springer.
