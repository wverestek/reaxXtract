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
    else:
        raise ValueError(f"File reader: Format {informat} for {infile} not yet supported!")


#########################
# read reaxff bond file #
#########################
def read_reax(infile: str) -> list[list[int,], list[nx.Graph,]]:
    # Initilize variables
    idx = -1
    ts = []
    pnum = []
    nxg = []

    # open file
    log.info(f"Reading reax file: {infile}")
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