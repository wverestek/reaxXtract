# reaXtract
reaXtract is a python tool to extract reactions from molecular dynamics
simulations with the reax force field.

At the moment this is just some rough code I put together and far from being 
finished or fully capable of an automatic analysis, reaction rates etc.

What reaXtract does:
- extracting reactions from reaxff bond files by identifying changed connectivity
- print reacting atom IDs to the console
- For each found Reaction (e.g. R0) two pictures are created in subfolder 'reaxff.bond.file'.dir
  
What reaXtract does not (maybe yet):
- take into account positions
- take into account unbonded interactions (catalysts, etc.)
- generate SMILES, SMART etc. output

## Citing

no publication yet. For the time being, please cite Wolfgang Verestek with the github repo

## Installation

no installation needed

## Prerequisites

Python modules:
- sys, os, random
- networkX

## Usage

python3 reaXtract.py 'reaxff.bonds.file'

stepwidth for frames can be adjusted at the beginning of the script:
- startstep: index of step to start with
- checkstep: compare to index+checkstep
- framestep: move index framestep every time.
e.g. startstep 0, checkstep 1, framestep 5: compare frame 0 with 1, 5 with 6, 10 with 11 etc.

colorcoding for the plots can be spcefied for each atom type in the dictionary type2color.
For colors see https://matplotlib.org/3.1.1/gallery/color/named_colors.html

## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.
