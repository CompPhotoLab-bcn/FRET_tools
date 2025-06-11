## FRET tools: fret_exp.py & fret_theo.py
This Python3 scripts implement the calculation of FRET efficiencies from fluorescence spectra (donor quenching at different acceptor concentrations) as well as from theoretical trajectories of electronic couplings. This software was developed in the framework of the ALLODD ITN project: https://www.allodd-itn.eu to screen collections of fluorescent small molecules to detect binding and analyze FRET data in protein-ligand complexes, FRET allowing to recover spatial information on the binding mode/site.

`fret_exp.py` estimates FRET efficiency from the quenching of the protein fluorescence (donor) at a given wavelength (with minimal acceptor emission) indicated in the input. It provides a plot in PDF format showing the fluorescence spectra at different ligand (acceptor) concentrations given as input. It also fits the fluorescence quenching data to the linear equation log(F0/F - 1) = log(Kb) + n log([L]), allowing to derive the association and dissociation constants (Kb=1/Kd) and ligand molecularity (n) from the fluorescence in the absence (F0) and presence (F) of the ligand at different concentrations [L]. Finally, it uses the association constant to estimate the percentage of free protein versus complex, in order to provide corrected FRET efficiencies taking into account the actual fraction of protein-ligand complex formed. All output results are printed in the output and saved to the fret_exp_out.csv file.

`fret_theo.py` estimates theoretical FRET efficiencies from electronic coupling trajectories obtained from a Molecular Dynamics trajectory of the protein-ligand complex. It reads an input file with several datasets of electronic coupling trajectories (e.g. those obtained in different binding sites for a given ligand, or for different MD replicas). The program computes the FRET efficiency in an intermediate coupling regime, where coupling fluctuations along the trajectory can have a similar timescale compared to that of the donor excited-state decay. It then also computes more approximate FRET efficiencies in the dynamic (coupling fluctuations fast compared with the donor excited-state decay) and static (coupling fluctuations slow compared with the donor excited-state decay) limits. The timescale of the donor excited-state decay is estimated from the donor fluorescence lifetime (default value 5.78 ns, measured for human serum albumin). The code also needs as input the donor-acceptor spectral overlap value (between normalized donor emission and normalized acceptor absoroption) needed in Förster FRET golden rule expression, as well as the timestep between coupling values in the input trajectories. 

Electronic coupling trajectories can be computed using Förster point dipole approximation, or from more rigorous calculations using QM-derived transition densities or transition charges. 

## References
Information on the theoretical background can be found in these publications:
L. Cupellini, M. Corbella, B. Mennucci and C. Curutchet, Electronic energy transfer in biomacromolecules, WIREs Comput. Mol. Sci. 2019, 9(2), e1392. https://doi.org/10.1002/wcms.1392

D. Gonzalo, L. Cupellini and C. Curutchet, On the breakdown of Förster energy transfer theory due to solvent effects: atomistic simulations unveil distance-dependent dielectric screening in calmodulin, Chem. Sci. 2025, 16(8), 3693–3704. https://doi.org/10.1039/D4SC07679F
 
## Disclaimer and Copyright
The terms for using, copying, modifying, or distributing this code are specified in the file `LICENSE`

## Authors
**Carles Curutchet**

**Özge Ergün**

Departament de Farmàcia i Tecnologia Farmacèutica, i Fisicoquímica & Institut de Química Teòrica i Computacional 

Facultat de Farmàcia i Ciències de l'Alimentació, Universitat de Barcelona

Av. Joan XXIII 27-31, 08028 Barcelona, Spain

## System requirements
These codes are implemented with python3, and require the following python libraries.

`fret_exp.py`:  `numpy`, `matplotlib` and `pandas`

`fret_theo.py`:  `numpy`

## Tests
The folder `src` contains the `fret_exp.py` and `fret_theo.py` codes.
The folder `test` contains a tests with inputs/outputs for each code:

Test for FRET_EXP can be executed running
```
cd test/test_exp
../../src/fret_exp.py naproxen_280nm.in naproxen_280nm.pdf --x_title "Emission wavelength (nm)" --y_title "Intensity" --plot_title "Naproxen (excitation at 280nm)" --width 7.0 --height 3.5 --wavelength 315 --show_plot True
```

Test for FRET_EXP can be executed running
```
cd test/test_theo
../../src/fret_theo.py --lifetime 5.78 --overlap 0.9773 --step 49 --coup coup_HSA.in --out HSA
```

## Usage 
Help on the usage of each code can be obtained with these commands
```
src/fret_exp.py --help
src/fret_theo.py --help
```

FRET_EXP options
```
usage: fret_exp.py [-h] [-x X_TITLE] [-y Y_TITLE] [-p PLOT_TITLE] [-W WIDTH] [-H HEIGHT] [-l WAVELENGTH] [-s SHOW_PLOT] dataset_info output_file

Compute FRET efficiency from fluorescence quenching data

positional arguments:
  dataset_info          Path to input file indicating flu datasets to be read
  output_file           Path to output PDF file with fluorescence plot

options:
  -h, --help            show this help message and exit
  -x X_TITLE, --x_title X_TITLE
                        Title for X-axis (default Emission wavelength (nm))
  -y Y_TITLE, --y_title Y_TITLE
                        Title for Y-axis (default Intensity)
  -p PLOT_TITLE, --plot_title PLOT_TITLE
                        Title for the plot (default Fluorescence spectra)
  -W WIDTH, --width WIDTH
                        Width of the plot in inches (default 14.0)
  -H HEIGHT, --height HEIGHT
                        Height of the plot in inches (default 7.0)
  -l WAVELENGTH, --wavelength WAVELENGTH
                        Wavelength at which intensities/FRET efficiencies will be computed
  -s SHOW_PLOT, --show_plot SHOW_PLOT
                        Whether to display the plot (detault True)
```
The dataset_info input is a plain text file indicating the names of the files with fluorescence spectra to be read, and the donor (protein) and acceptor (ligand) concentrations for each one in  μM units. The program needs one spectra with 0 ligand (acceptor) concentration as reference. For example, the `naproxen_280nm.in` file in the test provided indicates a reference HSA (human serum albumin) spectra at 5 μM and HSA-naproxen spectra at 10-100 μM ligand concentrations:
```
                        [HSA] [1N]
HSA5_280nm.txt            5     0
HSA5+NAP10_280nm.txt      5    10
HSA5+NAP50_280nm.txt      5    50
HSA5+NAP100_280nm.txt     5   100
```
Each file `HSA5_280nm.txt` etc contains the raw data of a fluorescence spectra with X Y format indicating Emission wavelength and Intensity on the first two columns.

FRET_THEO options
```
usage: fret_theo.py [-h] [--lifetime LIFETIME] --overlap J [--step TIMESTEP] [--coup COUP] [--out OUT]

Compute FRET efficiencies from multiple electronic coupling trajectories

options:
  -h, --help           show this help message and exit
  --lifetime LIFETIME  Donor fluorescence lifetime in ns (default 5.78 ns)
  --overlap J          Spectral overlap value in eV-1
  --step TIMESTEP      Time step in ps between coupling values in datasets (default 49 ps)
  --coup COUP          Electronic couplings file (in cm -1). Can contain mutliple columns/trajectories, with title in first line
  --out OUT            Prefix for output files with rate (ns-1) and efficiency distributions at each time t and time window
```
