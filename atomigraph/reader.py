# functions for handling data import and export
import os.path
import gzip
import glob
import re
from typing import List, Tuple, Union

import networkx as nx
from .logger import log
# from networkx.classes import Graph
# from .utils import ON2ELEM, ON2HEX, ELEM2HEX, DEFAULT_COLOR


def _natural_key(path: str):
    """
    Natural sort key for file paths similar to `ls -v`.
    Splits the basename into alternating non-digit and digit parts and
    returns a tuple where digit parts are ints (so 8 < 10).
    Works independent of differing name prefixes.
    """
    name = os.path.basename(path)
    parts = re.findall(r'\d+|\D+', name)
    key = []
    for p in parts:
        if p.isdigit():
            # preserve numeric ordering, handle big integers
            key.append(int(p))
        else:
            # case-insensitive text ordering
            key.append(p.lower())
    return tuple(key)


def _expand_infiles(infile: Union[str, List[str]]) -> List[str]:
    """
    Accepts a filename, a list of filenames or glob pattern(s),
    returns a version-aware sorted, validated list of existing file paths.
    Raises FileNotFoundError if no files found.
    """
    if infile is None or infile == "":
        return []

    items: List[str] = []
    if isinstance(infile, (list, tuple)):
        raw_items = list(infile)
    else:
        raw_items = [infile]

    for it in raw_items:
        if any(ch in it for ch in ("*", "?", "[")):
            matches = glob.glob(it)
            if matches:
                items.extend(matches)
        else:
            items.append(it)

    # validate existence and remove duplicates (preserve first occurrence)
    seen = set()
    unique_files: List[str] = []
    for p in items:
        if p in seen:
            continue
        if not os.path.exists(p):
            raise FileNotFoundError(f"Input file not found: {p}")
        seen.add(p)
        unique_files.append(p)

    if not unique_files:
        raise FileNotFoundError(f"No input files found for pattern(s): {infile}")

    # perform version-aware (natural) sort across all resolved files
    unique_files.sort(key=_natural_key)
    return unique_files

##########################
# read bond file wrapper #
##########################
def read_bonds(infile: Union[str, List[str]] = "", informat: str = "reaxff") -> list[list[int,], list[nx.Graph,]]:
    """
    Read one or multiple bond files.
    - infile: single filename, list of filenames, or glob pattern(s)
    - returns combined (ts, nxg) from all matched files (in version-aware sorted order)
    """
    files = _expand_infiles(infile)
    if informat.lower() == "reaxff":
        all_ts: List[int] = []
        all_nxg: List[nx.Graph] = []
        for f in files:
            ts, nxg = read_reax(f)
            all_ts.extend(ts)
            all_nxg.extend(nxg)
        return all_ts, all_nxg
    elif informat.lower() == "lammps_data":
        all_ts: List[int] = []
        all_nxg: List[nx.Graph] = []
        for f in files:
            ts, nxg = read_lammps_data(f)
            all_ts.extend(ts)
            all_nxg.extend(nxg)
        return all_ts, all_nxg
    else:
        raise ValueError(f"File reader: Format {informat} for {infile} not yet supported!")


#########################
# read reaxff bond file #
#########################
def read_reax(infile: str) -> list[int, nx.Graph]:
    # Initilize variables
    idx = -1
    ts = []
    pnum = []
    nxg = []

    # open file
    log.info(f"Reading reax file: {repr(infile)}")
    opener = gzip.open if infile.endswith(".gz") else open
    mode = "rt" if infile.endswith(".gz") else "r"
    with opener(infile, mode) as f:
        # iterate file line by line (more efficient than repeated readline calls)
        for line in f:
            if not line:
                break
            varline = line.strip().split()
            if line.startswith("# Timestep"):
                # Timestep
                idx = idx + 1
                nxg.append(nx.Graph())         # array of networkx graphs
                ts.append(int(varline[-1]))    # array of timesteps
                log.info(f"Reading {infile}\tFrame: {idx}\tTimestep: {ts[idx]}")
                continue
            elif line.startswith("# Number of particles"):
                # Number of particles
                pnum.append(int(varline[-1]))
                continue
            elif line.startswith("#") or line.startswith("\n") or len(varline) == 0:
                # other header lines or empty line
                continue
            elif varline[0].isdigit():
                # lines with atom/bond info
                aidx = 0  # atom index
                atomID = [0] * pnum[idx]    # atom ID
                atomType = [0] * pnum[idx]  # atom Type
                abo = [0.0] * pnum[idx]     # atom bond order
                nlp = [0.0] * pnum[idx]     # non-linarized potential
                q = [0.0] * pnum[idx]       # atom charge
                mol = [0] * pnum[idx]       # molecule ID
                bonds = []                  # bond list

                # process current line and subsequent atom lines using the file iterator
                cur_varline = varline
                while True:
                    log.log(5, f"line: {' '.join(cur_varline)}")
                    # atom info
                    atomID[aidx] = int(cur_varline[0])      # atom ID
                    atomType[aidx] = int(cur_varline[1])    # atom Type
                    abo[aidx] = float(cur_varline[-3])      # atom bond order
                    nlp[aidx] = float(cur_varline[-2])      # non-linarized potential
                    q[aidx] = float(cur_varline[-1])        # atom charge

                    # bond info
                    nb = int(cur_varline[2])
                    mol[aidx] = int(cur_varline[3 + nb])
                    for tmp in range(nb):
                        bonds.append((atomID[aidx], int(cur_varline[3 + tmp]), float(cur_varline[3 + nb + 1 + tmp])))

                    # advance to next line from the file iterator
                    next_line = next(f, None)
                    if next_line is None:
                        # EOF -> finish timestep processing
                        break

                    next_varline = next_line.strip().split()
                    # if header or empty line -> end of atom block for this timestep
                    if next_line.startswith("#") or next_line.startswith("\n") or len(next_varline) == 0:
                        # we've consumed the header/blank line; the outer for-loop will continue after this point
                        break
                    else:
                        # continue with next atom line
                        aidx += 1
                        cur_varline = next_varline

                # End of timestep, fill Graph with atoms and bonds for this timestep
                tmp = [(a, {"type": b}) for a, b in zip(atomID, atomType)]
                nxg[idx].add_nodes_from(tmp)

                # edges = bonds
                tmp = [(a, b, {"bo": c}) for a, b, c in bonds]
                nxg[idx].add_edges_from(tmp)

        # implicit f.close() via context manager
    return [ts, nxg]


#########################
# read LAMMPS data file #
#########################
def read_lammps_data(infile: str) -> list[list[int,], list[nx.Graph,]]:
    ts = 0                                                      # in case there is no timestep in data file
    opener = gzip.open if infile.endswith(".gz") else open      # unlikely to be gzipped, but ....
    mode = "rt" if infile.endswith(".gz") else "r"
    with opener(infile, mode) as f:
        print(f"Reading LAMMPS data file {repr(infile)}")
        line = f.readline()
        while line:
            #print("process: ",line.strip())
            # header section
            if "timestep =" in line:
                ts = int(line.strip().split()[-1])              # should be last bosition in standard data file
            elif line.startswith("#") or len(line.strip())==0:
                pass
            elif "atoms" in line:
                natoms = int(line.strip().split()[0])
            elif "atom types" in line:
                natomtypes = int(line.strip().split()[0])
            elif "bonds" in line:
                nbonds = int(line.strip().split()[0])
            elif "bond types" in line:
                nbondtypes = int(line.strip().split()[0])
            elif "angles" in line:
                nangles = int(line.strip().split()[0])
            elif "angle types" in line:
                nangletypes = int(line.strip().split()[0])
            elif "dihedrals" in line:
                ndihedrals = int(line.strip().split()[0])
            elif "dihedral types" in line:
                ndihedraltypes = int(line.strip().split()[0])
            elif "impropers" in line:
                nimpropers = int(line.strip().split()[0])
            elif "improper types" in line:
                nimpropertypes = int(line.strip().split()[0])
            elif any(x in line for x in ["xlo", "ylo", "zlo","avec","bvec","cvec","xy","xz","yz"]):
                pass
                
            # Masses
            elif "Masses" in line:
                _ = f.readline() # skip one empty line
                idx2mass = dict()
                # read natoms masses
                for idx in range(natomtypes):
                    words = f.readline().strip().split()
                    idx2mass[int(words[0])] = float(words[1])
                _ = f.readline() # skip one empty line
                
            # Coeffs
            elif any(line.strip().startswith(x) for x in ["Pair Coeffs"]):
                #print("Skiping: Pair Coeffs")
                for idx in range(natomtypes+2):             # skip Coeffs
                    _ = f.readline()
            elif any(line.strip().startswith(x) for x in ["Bond Coeffs"]):
                #print("Skipping: Bond Coeffs")
                for idx in range(nbondtypes+2):             # skip Coeffs
                    _ = f.readline()
            elif any(line.strip().startswith(x) for x in ["Angle Coeffs","BondBond Coeffs","BondAngle Coeffs"]):
                #print("Skipping: Angle Coeffs")
                for idx in range(nangletypes+2):            # skip Coeffs
                    _ = f.readline()
            elif any(line.strip().startswith(x) for x in ["Dihedral Coeffs","AngleAngleTorsion Coeffs","EndBondTorsion Coeffs",
                                        "MiddleBondTorsion Coeffs","BondBond13 Coeffs","AngleTorsion Coeffs"]):
                #print("Skipping: Dihedral Coeffs")
                for idx in range(ndihedraltypes+2):         # skip Coeffs
                    _ = f.readline()
            elif any(line.strip().startswith(x) for x in ["Improper Coeffs","AngleAngle Coeffs"]): 
                #print("Skipping: Improper Coeffs")
                for idx in range(nimpropertypes+2):         # skip Coeffs
                    _ = f.readline()
                
            # Atoms, bonds, Angles, Dihedrals, Impropers
            elif "Atoms" in line:
                _ = f.readline() # skip one empty line
                if "full" in line:
                    atomID = [[] for _ in range(natoms)]
                    mol = [[] for _ in range(natoms)]
                    atomType = [[] for _ in range(natoms)]
                    q = [[] for _ in range(natoms)]
                    pos = [[None,None,None] for _ in range(natoms)]
                    # skip the rest
                    for idx in range(natoms):
                        line = f.readline()
                        words = line.strip().split()
                        atomID[idx],mol[idx],atomType[idx] = [int(i) for i in words[0:3]]
                        q[idx],pos[idx][0],pos[idx][1],pos[idx][2] = [float(i) for i in words[3:7]]
                else:
                    #log.error("")
                    print("ERROR: unknown atom style")
                    sys.exit(0) 
                _ = f.readline() # skip one empty line
            elif "Velocities" in line:
                for idx in range(natoms+2): 
                    _ = f.readline() # skip Velocities
            elif "Bonds" in line:
                _ = f.readline() # skip one empty line
                bondID = [[] for _ in range(nbonds)]
                bondType = [[] for _ in range(nbonds)]    
                bondPair = [[None,None] for _ in range(nbonds)]
                for idx in range(nbonds):
                    line = f.readline()
                    words = line.strip().split()
                    bondID[idx],bondType[idx],bondPair[idx][0],bondPair[idx][1] = [int(i) for i in words]
                _ = f.readline() # skip one empty line
            elif "Angles" in line:
                for idx in range(nangles+2): 
                    _ = f.readline() # skip Angles
            elif "Dihedrals" in line:
                for idx in range(ndihedrals+2):
                    _ = f.readline() # skip Dihedrals
            elif "Impropers" in line:
                for idx in range(nimpropers+2): 
                    _ = f.readline() # skip Impropers
            else:
                #log.warning
                print("You should not be here!")
                print("##### ##### #####")
                print(line)
                print("##### ##### #####")
                sys.exit(0)
            # next line
            line = f.readline()
    
    # end of file
    
    # fill NetworkX Graph
    G = nx.Graph()
    # atoms => nodes
    tmp = [(a, {"type": b, "mass":c, "q":d}) for a,b,c,d in zip(atomID, atomType, [idx2mass[i] for i in atomType], q)]
    G.add_nodes_from(tmp)
    # edges => bonds
    tmp = [(a,b,{"bid":c,"bt":d}) for [a, b],c,d in zip(bondPair,bondID,bondType)]
    G.add_edges_from(tmp)
    
    return [[ts],[G]]
    
