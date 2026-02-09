#!/usr/bin/env python3

import sys
import argparse
from core import ReaxXtract

def main():
    parser = argparse.ArgumentParser(
        description="Analyze LAMMPS reax/bond output for bond topology changes."
    )
    parser.add_argument("-i", "--input", required=True, help="Input file (reax/bond format)")
    parser.add_argument("-a", "--atom-map", default="", help="Atom type mapping, e.g. '1:6,2:1,3:8'")
    parser.add_argument("-b", "--basename", default="", help="Output basename")
    parser.add_argument("--checkframe", type=int, default=1, help="Frame difference to check")
    parser.add_argument("--framestep", type=int, default=1, help="Frame step")
    parser.add_argument("--count-rings", action="store_true", help="Count ring structures")
    parser.add_argument("--ring-limits", default="3:10", help="Min:Max ring size")
    args = parser.parse_args()

    rxt = ReaxXtract(
        infile=args.input,
        basename=args.basename,
        atom_type_map=args.atom_map,
        checkframe=args.checkframe,
        framestep=args.framestep
    )

    try:
        rxt.read()
        rxt.find_rxns()
        rxt.plot_rxns()
        if args.count_rings:
            limits = tuple(map(int, args.ring_limits.split(":")))
            rxt.count_rings(ring_limits=limits)
    except Exception as e:
        print(f"Error: {e}")
        usage()
        sys.exit(2)

def usage():
    """Usage:
            python3 reaxXtract.py -i bonds.reaxff.dump -a=1:6,2:1,3:1,4:8,5:8,6:8,7:8,8:8 -r=3:10
            -i           input file in reax/bond format
            -a           atom type mapping (atom type : atomic number), e.g. 1:6 if atom type is carbon
            -b           basename for output, default: input filename without extension
            -r[=min:max] activates ring counting, 'min:max' atoms per ring, default 3:10
    """
    print(usage.__doc__)

if __name__ == "__main__":
    main()
    sys.exit(0)
