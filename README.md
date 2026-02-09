# reaxXtract
reaxXtract is a python tool to extract reactions from molecular dynamics
simulations with the reax force field (ReaxFF).

At the moment this is just some rough code I put together during my PhD thesis 
at IMWF University of Stuttgart and far from being professionally coded or 
fully capable of an automatic analysis, reaction rates etc... 
Provided as is, without any warranty, and with the hope that it will
be further developed in the future. Please use at your own risk, but feel free 
to use and modify it for your own purposes. Contributions welcome.

What reaxXtract does:
- extracting reactions from reaxff bond topology files by identifying changed connectivity
- store reactions in a pandas dataframe with time step, etc.
- remove reversed reactions within a certain time window (e.g. A+B -> C and C -> A+B)
- create a picture (before and after) for each reaction
  
What reaxXtract does not (yet?):
- take into account positions
- take into account non-explicit bond interactions (ions, hydrogen bridges, etc.)
- generate SMILES, SMART etc. output

## Citing

no publication yet. For the time being, please cite Wolfgang Verestek with the github repo

## Prerequisites

Python modules:
- sys, os, random
- numpy
- pandas
- matplotlib
- networkX
 
## Installation

add this folder to your PYTHONPATH.

## Usage

```
python3
import reaxXtract as rxt 
rxns=rxt.ReaxXtract(infile="bonds.reaxff.dump",atom_type_map="1:C,2:H,3:H,4:O,5:O,6:O,7:O,8:O")
rxns.read()
rxns.find_rxns()
# rxns.remove_reversed(timewindow=10)
rxns.plot_rxns()
```

command line usage:
```
tbd
python3 reaxXtract.py -i bonds.reaxff.dump -t 1:C,2:H,3:H,4:O,5:O,6:O,7:O,8:O
```
stepwidth for frames can be adjusted at the beginning of the script:
- startstep: index of step to start with
- checkstep: compare to index+checkstep
- framestep: move index framestep every time.
- e.g. startstep 0, checkstep 1, framestep 5: compare frame 0 with 1, 5 with 6, 10 with 11 etc.

Colorcoding for the plots can be spcefied for each atom type in the dictionary type2color in "utils.py".
Predefined colors follow Jmol color coding, see https://jmol.sourceforge.net/jscolors/

## Contributing

Contributions welcome
