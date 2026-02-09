# functions for handling data import and export
import os.path, gzip
import networkx as nx
from .logger import log
#from networkx.classes import Graph
#from .utils import ON2ELEM, ON2HEX, ELEM2HEX, DEFAULT_COLOR

##########################
# read bond file wrapper #
##########################
def read_bonds(infile:str="", informat:str="reaxff") -> list[ list[int,], list[nx.Graph,] ]:
    if not os.path.isfile(infile):
        raise FileNotFoundError(f"Input file not found: {infile}")
    
    if informat.lower() == "reaxff":
        [ts, nxg] = read_reax(infile)
    else:
        raise ValueError(f"File reader: Format {informat} for {infile} not yet supported!")
    return [ts, nxg]



#########################
# read reaxff bond file #
#########################
def read_reax(infile:str) -> list[ list[int,], list[nx.Graph,] ]:
    # open file
    log.info(f"Reading reax file: {infile}")
    opener = gzip.open if infile.endswith(".gz") else open
    mode = "rt" if infile.endswith(".gz") else "r"
    with opener(infile, mode) as f:
    
        # Initilize variables
        idx = -1
        ts = []
        pnum = []
        nxg = []
    
        # read file line by line
        line = f.readline()
        while line:
            if not line:
                break
            else:
                varline = line.strip().split()
                if line.startswith("# Timestep"):
                    # Timestep
                    idx = idx + 1
                    nxg.append(nx.Graph())         # array of networkx graphs
                    ts.append(int(varline[-1]))    # array of timesteps
                    log.info(f"Reading {infile}\tFrame: {idx}\tTimestep: {ts[idx]}")
                    line = f.readline()
                    continue
                elif line.startswith("# Number of particles"):
                    # Number of particles
                    pnum.append(int(varline[-1]))
                    line = f.readline()
                    continue
                elif line.startswith("#") or line.startswith("\n") or len(varline) == 0:
                    # other header lines or empty line
                    line = f.readline()
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
                    # read all lines with atom/bond info until next "#", "\n" or empty
                    while line:
                        log.log(5,f"line: {line.strip()}")
                        if not line:
                            # empty line, break while loop
                            break
                               
                        # atom info
                        atomID[aidx] = int(varline[0])      # atom ID
                        atomType[aidx] = int(varline[1])    # atom Type
                        abo[aidx] = float(varline[-3])      # atom bond order
                        nlp[aidx] = float(varline[-2])      # non-linarized potential
                        q[aidx] = float(varline[-1])        # atom charge
                    
                        # bond info
                        nb = int(varline[2])
                        mol[aidx] = int(varline[3 + nb])
                        for tmp in range(nb):
                            bonds.append((atomID[aidx], int(varline[3 + tmp]), float(varline[3 + nb + 1 + tmp])))
                    
                        # read next line
                        line = f.readline()
                        if line.startswith("#") or line.startswith("\n") or len(varline) == 0:
                            # break while loop, fill graph and continue next timestep
                            break
                        else:                     
                            # continue while loop with next line
                            aidx = aidx + 1
                            varline = line.strip().split()
                    
                    # End of timestep, fill Graph with atoms and bonds for this timestep
                    # nodes = atoms
                    #tmp = [(a, {"type": b, "abo": c, "nlp": d, "q": e,
                    #            "element":ON2ELEM.get(TYPE2ON.get(b,0),b),
                    #            "color":ON2HEX.get(TYPE2ON.get(b,0),DEFAULT_COLOR),
                    #            "mol": f})
                    #       for a, b, c, d, e, f in zip(atomID, atomType, abo, nlp, q, mol)]
                    tmp = [(a, {"type": b}) for a, b in zip(atomID, atomType)]
                    nxg[idx].add_nodes_from(tmp)    
                
                    # edges = bonds
                    tmp = [(a, b, {"bo": c}) for a, b, c in bonds]
                    nxg[idx].add_edges_from(tmp)
                
            # continue with next section "# Timestep" etc. or return
        f.close()
    return [ts, nxg]