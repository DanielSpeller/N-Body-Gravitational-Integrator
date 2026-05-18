#!/usr/bin/env python3
"""
N-Body Gravitational Integrator - Analysis Script

Reads the trajectory CSV and produces:
  1. A 3D trajectory plot (one trace per body)
  2. A 2D energy drift plot (relative energy error vs time)

Usage:
    python analysis/plot.py output/trajectory.csv [--save]

    --save  Write PNG files instead of showing interactive windows.
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def load_csv(path: str):
    return np.genfromtxt(path, delimiter=",", names=True)


def get_body_names(data) -> list:
    names = []
    for col in data.dtype.names:
        if col.endswith("_x"):
            names.append(col[:-2])
    return names


def plot_trajectories(data, body_names: list, save: bool = False):
    fig = plt.figure(figsize=(10, 8))
    ax  = fig.add_subplot(111, projection="3d")

    colors = plt.cm.tab10.colors

    for i, name in enumerate(body_names):
        x = data[f"{name}_x"]
        y = data[f"{name}_y"]
        z = data[f"{name}_z"]
        color = colors[i % len(colors)]
        ax.plot(x, y, z, label=name, color=color, linewidth=0.8, alpha=0.85)
        ax.scatter(x[0],  y[0],  z[0],  color=color, marker="o", s=60,  zorder=5)
        ax.scatter(x[-1], y[-1], z[-1], color=color, marker="*", s=120, zorder=5)

    ax.set_xlabel("x (AU)")
    ax.set_ylabel("y (AU)")
    ax.set_zlabel("z (AU)")
    ax.set_title("N-Body Trajectories")
    ax.legend(loc="upper left")
    plt.tight_layout()

    if save:
        plt.savefig("trajectories.png", dpi=150)
        print("Saved trajectories.png")
    else:
        plt.show()


def plot_energy_drift(data, save: bool = False):
    t      = data["time"]
    energy = data["energy"]
    E0     = energy[0]

    if E0 == 0.0:
        drift  = energy - E0
        ylabel = "Absolute energy drift  dE"
    else:
        drift  = (energy - E0) / abs(E0)
        ylabel = "Relative energy drift  (E - E0) / |E0|"

    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(t, drift, color="steelblue", linewidth=0.8)
    ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel(ylabel)
    ax.set_title("Energy Conservation - Velocity-Verlet Integrator")
    ax.grid(True, alpha=0.4)

    max_idx = int(np.argmax(np.abs(drift)))
    ax.annotate(
        f"max |dE/E0| = {drift[max_idx]:.2e}",
        xy=(t[max_idx], drift[max_idx]),
        xytext=(t[max_idx], drift[max_idx] * 1.5 if drift[max_idx] != 0 else 1e-9),
        fontsize=9,
        arrowprops=dict(arrowstyle="->", color="red"),
        color="red",
    )

    plt.tight_layout()

    if save:
        plt.savefig("energy_drift.png", dpi=150)
        print("Saved energy_drift.png")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description="Visualise N-body simulation output.")
    parser.add_argument("csv", help="Path to trajectory CSV file")
    parser.add_argument("--save", action="store_true",
                        help="Save plots as PNG files instead of displaying")
    args = parser.parse_args()

    data       = load_csv(args.csv)
    body_names = get_body_names(data)

    print(f"Loaded {len(data)} time steps, bodies: {body_names}")

    plot_trajectories(data, body_names, save=args.save)
    plot_energy_drift(data, save=args.save)


if __name__ == "__main__":
    main()
