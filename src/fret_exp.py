#!/usr/bin/env python
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Version
VERS = "1.0"

# Function to find the closest value to the target in a list
def find_closest_value(target, data):
    return min(data, key=lambda x: abs(x - target))

# Custom function to convert string argument to boolean
def str2bool(v):
    return v.lower() in ('true', 'yes', '1', 'y')

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
   print("                                   |     F R E T - E X P vers %-6.5s     |    " % VERS)
   print("                                   ╞ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ═ ╡    ")
   print("                                   | Carles Curutchet                    |    ")
   print("                                   | Ozge Ergun                          |    ")
   print("                                   | www.ub.edu/cplab                    |    ")
   print("                                   | www.allodd-itn.eu                   |    ")
   print("                                   ╰-------------------------------------╯    ")

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Compute FRET efficiency from fluorescence quenching data')
parser.add_argument('dataset_info', type=str, help='Path to input file indicating flu datasets to be read')
parser.add_argument('output_file', type=str, help='Path to output PDF file with fluorescence plot')
parser.add_argument('-x','--x_title', type=str, default='Emission wavelength (nm)', help='Title for X-axis (default Emission wavelength (nm))')
parser.add_argument('-y','--y_title', type=str, default='Intensity', help='Title for Y-axis (default Intensity)')
parser.add_argument('-p','--plot_title', type=str, default='Fluorescence spectra', help='Title for the plot (default Fluorescence spectra)')
parser.add_argument('-W','--width', type=float, default=14.0, help='Width of the plot in inches (default 14.0)')
parser.add_argument('-H','--height', type=float, default=7.0, help='Height of the plot in inches (default 7.0)')
parser.add_argument('-l','--wavelength', type=float, help='Wavelength at which intensities/FRET efficiencies will be computed')
parser.add_argument('-s','--show_plot', type=str, default='True', help='Whether to display the plot (detault True)')
args = parser.parse_args()

# Convert string argument to boolean
args.show_plot = str2bool(args.show_plot)

welcome()

# Read the dataset information from the text file
dataset_info = []
p_concentrations = set()  # To store unique p_concentration values
l_concentrations = set()  # To store unique l_concentration values
efficiency_data = []  # To store efficiency data for non-zero l_concentration datasets
F0_intensity = None

with open(args.dataset_info, 'r') as file:
    next(file)  # Skip the first line
    for line in file:
        filename, p_concentration, l_concentration = line.split()
        p_concentration = float(p_concentration)
        l_concentration = float(l_concentration)

        # Add to dataset_info for plotting (even if p_concentration == 0.0)
        dataset_info.append((filename, p_concentration, l_concentration))
        p_concentrations.add(p_concentration)
        l_concentrations.add(l_concentration)

# Check if there is exactly one dataset with 0 l_concentration and the rest have positive l_concentration
zero_l_concentration_datasets = sum(1 for l_concentration in l_concentrations if l_concentration == 0)
positive_l_concentration_datasets = len(l_concentrations) - zero_l_concentration_datasets

if zero_l_concentration_datasets != 1 or positive_l_concentration_datasets < 1:
    raise ValueError("There should be exactly one dataset with 0 l_concentration and at least one dataset with positive l_concentration.")

# Plotting setup
plt.figure(figsize=(args.width, args.height))

# Process each dataset for plotting
for info in dataset_info:
    filename, p_concentration, l_concentration = info

    # Read the dataset file
    x = []
    y = []
    with open(filename, 'r') as dataset_file:
        for line in dataset_file:
            x_val, y_val = line.split()
            x.append(float(x_val))
            y.append(float(y_val))

    # Plot the dataset with legend
    p_label = f'[P] = {p_concentration} μM'
    l_label = f'[L] = {l_concentration} μM'
    label = f'{p_label}, {l_label}'
    plt.plot(x, y, label=label)

    # Find and print the intensity value at the given wavelength
    if args.wavelength is not None:
        closest_wavelength = find_closest_value(args.wavelength, x)
        index = x.index(closest_wavelength)
        intensity_at_wavelength = y[index]

        # Handle calculations only for datasets with non-zero p_concentration
        if p_concentration > 0:
            # Calculate and store efficiency data for datasets with non-zero l_concentration
            if l_concentration == 0:
                F0_intensity = intensity_at_wavelength
                print("\nReference spectra (Protein donor only):")
                print(f"Ligand concentration (acceptor): {l_concentration} μM")
                print(f"Protein concentration (donor): {p_concentration} μM")
                print(f"Wavelength: {args.wavelength} nm")
                print(f"Reference intensity at {args.wavelength} nm {intensity_at_wavelength}")

                efficiency_data.append((l_concentration, p_concentration, args.wavelength, intensity_at_wavelength, 0))  # Efficiency is 0 for 0 l_concentration
            elif l_concentration > 0:
                efficiency = 1.0 - (intensity_at_wavelength / F0_intensity)
                efficiency_data.append((l_concentration, p_concentration, args.wavelength, intensity_at_wavelength, efficiency))

# Add labels and legend
plt.xlabel(args.x_title)
plt.ylabel(args.y_title)
plt.title(args.plot_title)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Adjust layout and spacing
plt.tight_layout()

# Save the plot as a PDF file
plt.savefig(args.output_file, format='pdf', bbox_inches='tight')

# Write efficiency data to the output file
if efficiency_data:
    df = pd.DataFrame(efficiency_data, columns=['l_concentration', 'p_concentration', 'wavelength', 'intensity_at_wavelength', 'efficiency'])
    df_nonzero_l_concentration = df[df['l_concentration'] > 0].copy()  # Filter out the data for 0 l_concentration

    # Calculate values
    x_values = np.log10(df_nonzero_l_concentration['l_concentration'])
    y_values = np.log10(F0_intensity / df_nonzero_l_concentration['intensity_at_wavelength'] - 1)

    # Perform linear regression for non-zero l_concentration data only
    slope, intercept = np.polyfit(x_values, y_values, 1)

    # Calculate correlation coefficient (r)
    correlation_matrix = np.corrcoef(x_values, y_values)
    r = correlation_matrix[0, 1]

    print("\nCalculation of binding constant: Linear regression data:")
    for x_val, y_val in zip(x_values, y_values):
        print(f"log([L]): {x_val:.4f}, log(F0/F - 1): {y_val:.4f}")

    # Calculate n and Kb
    n = slope
    Kb = 10 ** intercept
    Kd = 1 / Kb

    print("\nEquation to fit:")
    print(f"log(F0/F - 1) = log(Kb) + n log([L])")

    print("\nLinear regression results:")
    print(f"Slope (n, binding molecularity): {n:.4f}")
    print(f"Intercept (log(Kb)): {intercept:.4f}")
    print(f"Association constant Kb: {Kb:.4f} μM⁻¹")
    print(f"Dissociation constant Kd: {Kd:.4f} μM")
    print(f"Correlation coefficient (r): {r:.4f}")

    # Calculate percent_free_prot using the formula 1 / (Kb * l_concentration + 1)
    percent_free_prot = 1 / (Kb * df_nonzero_l_concentration['l_concentration'] + 1)

    # Calculate percent_complex as 1 - percent_free_prot
    percent_complex = 1 - percent_free_prot

    # Calculate efficiency_corrected as efficiency / percent_complex
    efficiency_corrected = df_nonzero_l_concentration['efficiency'] / percent_complex

    # Add calculated values to the DataFrame using .loc
    df_nonzero_l_concentration.loc[:, 'percent_free_prot'] = percent_free_prot
    df_nonzero_l_concentration.loc[:, 'percent_complex'] = percent_complex
    df_nonzero_l_concentration.loc[:, 'efficiency_corrected'] = efficiency_corrected

    print("\nFRET efficiency uncorrected/corrected by percent of free protein vs complex:")
    print(df_nonzero_l_concentration)

    # Write efficiency data to the output file
    with open('fret_exp_out.csv', 'w') as output_table:
        df_nonzero_l_concentration.to_csv(output_table, index=False)

else:
    print("No datasets with non-zero ligand concentration found.")

# Show the plot if required
if args.show_plot:
    plt.show()

