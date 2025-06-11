#!/usr/bin/env python
import numpy as np
import argparse
import os

# Version
VERS = "1.0"

# Constants
HBAR = 6.582119569e-16  # Planck's hbar constant in eV*s
TWO_PI = 2*np.pi
EV2WN = 8065.54477

def welcome():
   print("                                                                              ")
   print("                                                                  ▓▒▒▒▒       ")
   print("                                                      ░▓▓▒▒▒    ▒▓▒▒▒▒░       ")
   print("                                ░▓▒▒▒▒    ▓▓▒▒░░░    ▓▒▒░░░░░   ░░▒▒▒  ▒      ")
   print("           ░▒░░░░░    ▒▒░░░░    ▒▒░░░░░   ░░▒░░░░░   ░░░▒░░░░░ ░░░░░     ▒    ")
   print("             ░░░░    ░░░░░░░   ░░░▒░░░░  ░░░░▒░░░░░ ░░░░░▒░░░░░▒░░░           ")
   print("              ░░░░░  ░░░▒░░░░░ ░░░░▒░░░░░░░░░ ▒▒░░▒░▒░░░  ▓▒▒▒▓▓▓▒▒           ")
   print("              ▒▒▒▒▒░░░░░ ▒▒▒▒▒▒░░░░ ▓▒▒▒▒▒░░░  ▒▒▒▒▓▓▒▒░                      ")
   print("           ░    ▓▓▓▓▓▒▒░   ▓▓▓▓▓▒▒                                            ")
   print("   ░     ____                      ____  _           _        _          _    ")
   print("       / ___|___  _ __ ___  _ __ |  _ \| |__   ___ | |_ ___ | |    __ _| |__  ")
   print("      | |   / _ \| '_ ` _ \| '_ \| |_) | '_ \ / _ \| __/ _ \| |   / _` | '_ \ ")
   print(" ░    | |__| (_) | | | | | | |_) |  __/| | | | (_) | || (_) | |__| (_| | |_) |")
   print("       \____\___/|_| |_| |_| .__/|_|   |_| |_|\___/ \__\___/|_____\__,_|_.__/ ")
   print("   ░                       |_|                                                ")
   print("    ░                                                        ▓▓▒▒▒▒▒          ")
   print("      ▒░░░ ▒                                        ░▓▒▒▒  ░▒▒▒▒░   ░         ")
   print("          ░░░░░                ▓▒▒     ▓▓▓▒▒▒    ▒▓▓▒░░░░░ ░░░░░     ▒ ▒      ")
   print("          ░░░░     ▒▒▒▒▒▒    ▒▒▒░░░    ▒▒▒░░░░   ░░░░▒░░   ░░  ░      ▒▒      ")
   print("          ░░░░    ░░░░░░░░  ░░▒░░░░   ░░░░░░░░  ░░░░░░ ░░░░░░░░               ")
   print("          ▒▒▒░░ ░░░░░░░  ░  ░░░▒▒░░░  ░░░░░░░░░░░░░░    ▒▒▒▒▒▒▓               ")
   print("           ░▓▓▒▒░░░ ░░░░░░ ░░░  ▒▒▒▒▒▒▒░░░ ▒▒▒▒▒▒▒▒░░                         ")
   print("                ▓░░  ▒░░░░░░░░     ▒▓▓▓░    ░▓▓▓▒                             ")
   print("                     ▓▒▒▒▒▒▒░                                                 ")
   print("                                   ╭-------------------------------------╮    ")
   print("                                   |     F R E T - T H E O  vers %-6.5s  |    " % VERS)
   print("                                   ╞ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ╡    ")
   print("                                   | Carles Curutchet                    |    ")
   print("                                   | Ozge Ergun                          |    ")
   print("                                   | www.ub.edu/cplab                    |    ")
   print("                                   | www.allodd-itn.eu                   |    ")
   print("                                   ╰-------------------------------------╯    ")

def compute_rate(V2, J):
    return (TWO_PI / HBAR) * (V2/EV2WN**2) * J * 1e-9

def compute_efficiency(lifetime, rate):
    return 1 / (1 + 1/(lifetime * rate))

def main():
    parser = argparse.ArgumentParser(description="Compute FRET efficiencies from multiple electronic coupling trajectories")
    parser.add_argument("--lifetime", type=float, default=5.78, help="Donor fluorescence lifetime in ns (default 5.78 ns)")
    parser.add_argument("--overlap", dest="J", type=float, required=True, help="Spectral overlap value in eV-1")
    parser.add_argument("--step", dest="timestep", type=float, default=49.0, help="Time step in ps between coupling values in datasets (default 49 ps)")
    parser.add_argument("--coup", default="coup.in", help="Electronic couplings file (in cm -1). Can contain mutliple columns/trajectories, with title in first line")
    parser.add_argument("--out", type=str, help="Prefix for output files with rate (ns-1) and efficiency distributions at each time t and time window")
    args = parser.parse_args()

    welcome()

    # Load the dataset from the specified couplings file
    data = np.loadtxt(args.coup, dtype=str)
    dataset_names = data[0]
    datasets = data[1:].astype(float)

    steps_per_lifetime = int( (args.lifetime*1e3) / args.timestep)

    for idx, name in enumerate(dataset_names):
        V2 = datasets[:, idx]**2
        rates_t = compute_rate(V2, args.J)
        eff_t = compute_efficiency(args.lifetime, rates_t)
        rate_win_values = []
       # Compute rate_win values by averaging over windows of size 'steps_per_lifetime'
        for i in range(0, len(rates_t) - steps_per_lifetime, 1):
            segment = rates_t[i:i+steps_per_lifetime]
            rate_win = np.mean(segment)
            rate_win_values.append(rate_win)

        eff_win = []

        for i in range(0, len(rate_win_values), 1):
            tmpeff = compute_efficiency(args.lifetime, rate_win_values[i])
            eff_win.append(tmpeff)

        efficiency = np.mean(eff_win)
        dynamic_rate = np.mean(rates_t)
        dynamic_efficiency = compute_efficiency(args.lifetime, dynamic_rate)
        static_efficiency = np.mean(eff_t)

        print(f"Results for {name}")
        print(f"   FRET efficiency (intermediate averaging regime) = {efficiency:.4e}")
        print(f"   Approximate results in the dynamic and static limits:")
        print(f"   FRET efficiency (dynamic averaging regime) = {dynamic_efficiency:.4e}")
        print(f"   FRET efficiency (static averaging regime) = {static_efficiency:.4e}")

        if args.out:
            prefix = f"{args.out}_{name}"
            np.savetxt(f"{prefix}_rates_win.txt", rate_win_values, fmt="%.4e")
            np.savetxt(f"{prefix}_eff_win.txt", eff_win, fmt="%.4e")
            np.savetxt(f"{prefix}_rates_t.txt", rates_t, fmt="%.4e")
            np.savetxt(f"{prefix}_eff_t.txt", eff_t, fmt="%.4e")

if __name__ == "__main__":
    main()

