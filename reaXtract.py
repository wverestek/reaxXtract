#!/usr/bin/env python3
import os.path
import random
import sys
import getopt
from getopt import GetoptError
from typing import TextIO

import matplotlib.pyplot as plt
import networkx as nx
from numpy.f2py.rules import typedef_need_dict


##########################
# convert string to dict # e.g. "1:1,2:6,3:8" -> {1:1,2:6,3:8}
##########################
def convert_str2dict(atom_type_map: str) -> dict:
    if len(atom_type_map) > 0 and ":" in atom_type_map:
        return {k: v for k, v in [[int(i) for i in tmp.split(":")] for tmp in atom_type_map.strip().split(",")]}
    else:
        return dict()

#############################################
# get k nearest neighbors of networkx Graph #
#############################################
def k_nearest_neighs(G:nx.Graph, start:set, k:int) -> set:
    neighs = start
    for l in range(k):
        neighs = set((nbr for n in neighs for nbr in G[n]))
    return neighs



class ReaXtract:
    ##############
    # initialize #
    ##############
    def __init__(self, infile:str="", intype:str="reaxff", basename:str="",
                 startstep:int=0, stopstep:int=sys.maxsize, checkframe:int=1, framestep:int=1,
                 hash_by:str="type", rxn_bond_cutoff:int=1, plot_bonds:int=5, seed:int=42,
                 atom_type_map:str="",
                 ring_counter=False,ring_limits=(3,10)):
        """
        A class to read reaxff bond information files  and extract changes in bond configuration

        infile : str
            A bond information file, that should be read. can also be set in class function read()
        intype : str
            file type of the file containing bond information. Default "reaxff"
        basename : str
            base name for output. If not set the input file name is used as base name.
        startstep : int
            MD time step to start evaluating the bond information. Default: 0
        stopstep : int
            MD time step to stop evaluating bond information. Default: sys.maxsize (a very high umber)
        checkframe : int
            Number of frames difference to check for changed bonds. Default: 1
        framestep : int
            Number of frames before the next evaluation is done. Default 1:
         e.g. If the bond file contains bond information every 10 steps from 0 to 100, the setting 0,100,1,1 will start
         at step 10 and compare it with 0, then 20 with 10 etc. A setting of 0,100,2,5 will start at frame 20, compare
         with 0, then continue with 70, compare to 50 , etc.
        hash_by : str
         hash for each reaction that is build upon 'element' or 'type'. Hashes allow to identify if a similar reaction
         has already occured before or not.

        """
        
        print("Initializing...")
        self.name:str = "reaXtract"
        self.infile:str = infile
        self.intype:str = intype.lower()
        if len(basename) > 0:
            self.basename:str = basename
        else:
            self.basename:str = infile

        self.startstep:int = int(startstep)
        self.stopstep:int = int(stopstep)
        self.startidx:int = 0
        self.stopidx:int = 0
        self.checkframe:int = int(checkframe)
        self.stepframe:int = int(framestep)

        self.ts:list = []
        self.nxg:list = []
        self.rxn_sets_before:list = [[]]
        self.rxn_sets_after:list = [[]]
        self.rxn_hashes_before:list = [[]]
        self.rxn_hashes_after:list = [[]]
        self.hash_by:str = hash_by
        self.rxn_elements_before = [[]]
        self.rxn_elements_after = [[]]
        self.rxn_id = []
        self.rxn_count = []

        self.ring_counter = ring_counter
        self.ring_limits = ring_limits
        self.ring_sets = []

        self.rxn_bond_cutoff:int = int(rxn_bond_cutoff)     # number of bonds distance that are considered one reaction, default: 0
        self.plot_bonds:int = int(plot_bonds)               # number of bonds distance to include for plots
        random.seed(a=int(seed))                            # reproducibility of pyplot plots
        self.default_color = "#E0E0E0"

        # atomic ordinal number to element dict
        self.on2elem: dict = {1: "H", 2: "He",
                              3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
                              11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar"}
        # atomic ordinal number to rgb hex color
        # antiquewhite FAEBD7 instead of white FFFFFF for Hydrogen
        self.on2rgb: dict = {1: "#FAEBD7", 2: "#D9FFFF", 3: "#CC80FF", 4: "#C2FF00", 5: "#FFB5B5", 6: "#909090",
                             7: "#3050F8", 8: "#FF0D0D", 9: "#90E050", 10: "#B3E3F5", 11: "#AB5CF2", 12: "#8AFF00",
                             13: "#8FAF00", 14: "#F0C8A0", 15: "#FF8000", 16: "#FFFF30", 17: "#1FF01F", 18: "#80D1E3"}

        # process some things
        # atom mapping from type to ordinal number -> user input as str
        self.type2on:dict = convert_str2dict(atom_type_map)
        self.elem2rgb:dict = {self.on2elem[i]:self.on2rgb[i] for i in range(1,len(self.on2elem)+1)}



    ##########################
    # read bond file wrapper #
    ##########################
    def read(self, infile:str="", intype:str="reaxff"):
        #print("read:",infile,self.infile)
        if len(infile)>0:
            self.infile = infile
        if len(self.basename) == 0 and len(self.infile)>0:
            self.basename = self.infile
        if intype.lower() == "reaxff":
            self.read_reax()
        else:
            print("File reader: intype", str(intype), "not yet supported!")

    #########################
    # read reaxff bond file #
    #########################
    def read_reax(self, infile:str=""):
        #print("read_reax:", infile,self.infile)
        if len(infile) > 0:
            self.infile = infile
        if len(self.basename) == 0 and len(self.infile)>0:
            self.basename = self.infile
            

        print("Reading reax file:", self.infile)
        f = open(self.infile, "r")
        idx = -1
        self.ts = []
        pnum = []
        self.nxg = []
        line = f.readline()
        while line:
            varline = line.strip().split()
            if line.startswith("# Timestep"):
                # Timestep
                idx = idx + 1
                self.nxg.append(nx.Graph())  # array of networkx graphs
                self.ts.append(int(varline[-1]))  # array of timesteps
                print("\rReading", self.infile, "\tFrame:", idx, "\tTimestep:", self.ts[idx], "...", end="")
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
                atomID = [0] * pnum[idx]
                atomType = [0] * pnum[idx]
                abo = [0.0] * pnum[idx]
                nlp = [0.0] * pnum[idx]
                q = [0.0] * pnum[idx]
                mol = [0] * pnum[idx]
                bonds = []
                # read all lines with atom/bond info until next "#", "\n" or empty
                while line:
                    # print(line)
                    # atom info
                    atomID[aidx] = int(varline[0])
                    atomType[aidx] = int(varline[1])
                    abo[aidx] = float(varline[-3])
                    nlp[aidx] = float(varline[-2])
                    q[aidx] = float(varline[-1])
                    # bond info
                    nb = int(varline[2])
                    mol[aidx] = int(varline[3 + nb])
                    for tmp in range(nb):
                        bonds.append((atomID[aidx], int(varline[3 + tmp]), float(varline[3 + nb + 1 + tmp])))

                    line = f.readline()
                    if line.startswith("#") or line.startswith("\n") or len(varline) == 0:
                        # break while loop, fill graph and continue next timestep
                        break
                    else:
                        # continue while loop with next line
                        aidx = aidx + 1
                        varline = line.strip().split()

                # End of timestep, fill Graph with atoms and bonds for this timestep
                # nodes
                tmp = [(a, {"type": b, "abo": c, "nlp": d, "q": e,
                            "element":self.on2elem.get(self.type2on.get(b,0),b),
                            "color":self.on2rgb.get(self.type2on.get(b,0),self.default_color)})
                       for a, b, c, d, e in zip(atomID, atomType, abo, nlp, q)]
                self.nxg[idx].add_nodes_from(tmp)
                # edges
                tmp = [(a, b, {"bo": c}) for a, b, c in bonds]
                self.nxg[idx].add_edges_from(tmp)

                # check for start and stop step# monotony in timesteps assumed
                if self.startstep == self.ts[-1]: self.startidx = idx
                if self.stopstep >= self.ts[-1]: self.stopidx = idx

                # continue with next section "# Timestep" etc. or return

            else:
                print("You shouldn't be here...")
        # finished reading file
        print(" finished")
        f.close()



    ##################
    # find reactions #
    ##################
    def find_rxn(self):
        print("Searching reactions...")
        f_rxn:TextIO = open(self.basename + "_rxnIDs.dat", "w")
        mystr = "# Timestep \t RxnID \t RxnCount \t FromIDs:ToIDs \t FromType:ToType \t FromElem:ToElem \t Rxn_hashes"
        print(mystr)
        f_rxn.write(mystr+"\n")

        start = self.startidx
        stop = self.stopidx
        cf = self.checkframe
        fs = self.stepframe

        self.rxn_sets_before = [[] for _ in range(len(self.ts))]
        self.rxn_sets_after =  [[] for _ in range(len(self.ts))]
        self.rxn_hashes_before = [[] for _ in range(len(self.ts))]
        self.rxn_hashes_after = [[] for _ in range(len(self.ts))]
        self.rxn_elements_before = [[] for _ in range(len(self.ts))]
        self.rxn_elements_after = [[] for _ in range(len(self.ts))]
        self.rxn_id = []
        self.rxn_count = []

        for idx in range(start + cf, stop+1, fs):
            # comparing graphs,
            # edges that are not in both graphs, list of tuples
            #print(idx-cf, idx)
            #print(nx.difference(self.nxg[idx - cf], self.nxg[idx]).edges())
            #print(nx.difference(self.nxg[idx], self.nxg[idx - cf]).edges())
            reacting_edges = set(nx.difference(rxt.nxg[idx -cf], rxt.nxg[idx]).edges()).union(
                             set(nx.difference(rxt.nxg[idx], rxt.nxg[idx - cf]).edges())       )
            # atom IDs that are involved with changed bond connectivity
            reacting_atoms = set(i for j in reacting_edges for i in j)                                       # set(int)

            # reactions found
            if len(reacting_atoms) > 0:
                # dicts for conversion
                node2element = nx.get_node_attributes(self.nxg[idx], name="element")
                node2type = nx.get_node_attributes(self.nxg[idx], name="type")

                # include reactions rxn_bond_cutoff away?
                print("0:",reacting_atoms)
                if self.rxn_bond_cutoff > 0:
                    reacting_atoms = k_nearest_neighs(self.nxg[idx-cf],reacting_atoms,self.rxn_bond_cutoff).union(
                                     k_nearest_neighs(self.nxg[idx]   ,reacting_atoms,self.rxn_bond_cutoff)       )  # combine sets
                print("1:",reacting_atoms)

                # group by connectivity for all reactions #
                # edges before and after reaction
                Gcombine = self.nxg[idx].subgraph(reacting_atoms).copy()                                   # nx.Graph
                before_edges = self.nxg[idx - cf].subgraph(reacting_atoms).edges()                         # nx.EdgeView([]) = list(tuple(int,int),), might be empty
                if len(before_edges) > 0:
                    Gcombine.add_edges_from(before_edges)
                # list of sets, one set per found reaction
                reaction_sets = list(nx.connected_components(Gcombine))                                    # list(set(int,),)
                del Gcombine


                # before and after sets for each individual reaction
                for rset in reaction_sets:
                    #print("rset:",rset)
                    # before
                    tmp = self.nxg[idx - cf].subgraph(rset)
                    connect1 = sorted(nx.connected_components(tmp))
                    hash1 = nx.weisfeiler_lehman_graph_hash(tmp, node_attr=self.hash_by)
                    self.rxn_sets_before[idx].append(connect1)
                    self.rxn_hashes_before[idx].append(hash1)
                    # after
                    tmp = self.nxg[idx].subgraph(rset)
                    connect2 = sorted(nx.connected_components(tmp))
                    hash2 = nx.weisfeiler_lehman_graph_hash(tmp, node_attr=self.hash_by)
                    self.rxn_sets_after[idx].append(connect2)
                    self.rxn_hashes_after[idx].append(hash2)

                # print output
                for i, (rxn_before, rxn_after) in enumerate(zip(self.rxn_sets_before[idx], self.rxn_sets_after[idx])):  # rxns: tuple[list[Any], list[Any]]
                    elem_before = [[] for _ in range(len(rxn_before))]
                    type_before = [[] for _ in range(len(rxn_before))]
                    # for each molecule before reaction
                    for j, atoms in enumerate(rxn_before):
                        # for each atom in molecule before reaction
                        for atom in list(atoms):
                            elem_before[j].append(node2element[atom])
                            type_before[j].append(node2type[atom])
                    elem_after = [[] for _ in range(len(rxn_after))]
                    type_after = [[] for _ in range(len(rxn_after))]
                    # for each molecule in after reaction
                    for j, atoms in enumerate(rxn_after):
                        # for each atom in molecule before reaction
                        for atom in list(atoms):
                            elem_after[j].append(node2element[atom])
                            type_after[j].append(node2type[atom])
                    self.rxn_elements_before[idx].append(elem_before)
                    self.rxn_elements_after[idx].append(elem_after)

                    rxn_hash = self.rxn_hashes_before[idx][i] + ':' + self.rxn_hashes_after[idx][i]
                    if rxn_hash not in self.rxn_id:
                        self.rxn_id.append(rxn_hash)
                        rxn_id = self.rxn_id.index(rxn_hash)
                        self.rxn_count.append(1)
                    else:
                        rxn_id = self.rxn_id.index(rxn_hash)
                        self.rxn_count[rxn_id] += 1
                    mystr = ";\t".join(tmp for tmp in
                            [str(self.ts[idx]),str(rxn_id),str(self.rxn_count[rxn_id]),
                             str(rxn_before)+":"+str(rxn_after),
                             str(type_before)+":"+str(type_after),
                             str(self.rxn_elements_before[idx][i])+":"+str(self.rxn_elements_after[idx][i]),
                             rxn_hash])
                    print(mystr)
                    f_rxn.write(mystr+"\n")


    # plot ractions with matplotlib
    def plot_rxn(self):
        print("plotting reactions...")
        start = self.startidx
        stop = self.stopidx
        cf = self.checkframe
        fs = self.stepframe

        # for each timestep
        for idx in range(start + cf, stop + 1, fs):

            # create output folder
            outfolder = self.basename + ".dir"
            if not os.path.exists(outfolder):
                os.makedirs(outfolder)

            node2elem = nx.get_node_attributes(self.nxg[idx], name="element")
            node2type = nx.get_node_attributes(self.nxg[idx], name="type")

            # for each reaction in timestep
            for i, rxn_before in enumerate(self.rxn_sets_before[idx]):
                rxn_atoms = set(i for j in rxn_before for i in j)
                edge_length = None #1.0/(len(rxn_atoms)**0.5)
                # get reacting edges
                G1 = self.nxg[idx-cf].subgraph(rxn_atoms)
                G2 = self.nxg[idx].subgraph(rxn_atoms)
                active_edges1 = nx.difference(G2, G1).edges()
                active_edges2 = nx.difference(G1, G2).edges()
                # reacting atoms and environment atoms
                print_atoms = k_nearest_neighs(self.nxg[idx], rxn_atoms, self.plot_bonds) # list of int
                Gprint = self.nxg[idx].subgraph(print_atoms).copy()
                # positions with spring model
                Gprint.add_edges_from(active_edges1)
                #pos = nx.forceatlas2_layout(Gprint,seed=random.randrange(0,5000))
                pos = nx.kamada_kawai_layout(Gprint)
                pos = nx.spring_layout(Gprint, iterations=200, pos=pos)

                del Gprint

                # before reaction
                plt.clf()
                Gprint = self.nxg[idx-cf].subgraph(print_atoms)
                color = [self.elem2rgb.get(v, self.default_color) for v in
                         nx.get_node_attributes(Gprint, name="element").values()]
                node_labels = {k:  str(node2elem[k]) + ":" + str(node2type[k]) + "\n" + str(k) for k in Gprint}
                bo = [tmp for tmp in list(nx.get_edge_attributes(Gprint, "bo").values())]
                pos = nx.spring_layout(Gprint, iterations=500, pos=pos)
                nx.draw_networkx_edges(Gprint, edgelist=active_edges2, alpha=0.6, width=5.0, edge_color="tab:red", pos=pos)
                nx.draw(Gprint, pos=pos, node_color=color, with_labels=True, labels=node_labels, font_size=6,
                        node_size=300, edge_color="black", width=bo)
                plt.tight_layout()
                f_out = os.path.join(outfolder, self.basename + "_" + str(self.ts[idx]) + "_Rxn" + str(i) + "_0.png")
                plt.savefig(f_out,dpi=200)

                # after reaction
                plt.clf()
                Gprint = self.nxg[idx].subgraph(print_atoms)
                color = [self.elem2rgb.get(v, self.default_color) for v in
                         nx.get_node_attributes(Gprint, name="element").values()]
                node_labels = {k:  str(node2elem[k]) + ":" + str(node2type[k]) + "\n" + str(k) for k in Gprint}
                bo = [tmp for tmp in list(nx.get_edge_attributes(Gprint, "bo").values())]
                pos = nx.spring_layout(Gprint, iterations=500, pos=pos)
                nx.draw_networkx_edges(Gprint, edgelist=active_edges1, alpha=0.6, width=5.0, edge_color="tab:red", pos=pos)
                nx.draw(Gprint, pos=pos, node_color=color, with_labels=True, labels=node_labels, font_size=6,
                        node_size=300, edge_color="black", width=bo)
                plt.tight_layout()
                f_out = os.path.join(outfolder, self.basename + "_" + str(self.ts[idx]) + "_Rxn" + str(i) + "_1.png")
                plt.savefig(f_out,dpi=200)
                # del reaction_sets



    # find and count rings #
    def count_rings(self, ring_limits:tuple=None):
        if self.ring_limits is None:
            self.ring_limits = (3, 10)
        if ring_limits != None:
            self.ring_limits = ring_limits

        print("Searching for ring structures with sizes:",self.ring_limits)
        start = self.startidx
        stop = self.stopidx
        cf = self.checkframe
        fs = self.stepframe
        ts = self.ts

        rl = self.ring_limits
        if len(rl) > 0:
            if len(rl) == 1 and rl[0] > 2:
                rl = (3, rl[0])
            elif len(rl) == 2:
                if rl[0] > rl[1]: rl = (rl[1], rl[0])
            else:
                print("ring count parameter ill defined")
                return 0

            # open file for output
            fring = open(self.basename + "_rc.dat", "w")
            a = "# timestep"
            b = " ".join(str(tmp) for tmp in list(range(rl[0], rl[1] + 1)))
            c = "total_ring_count"
            tmp = fring.write(" ".join([a, b, c, "\n"]))
            num_cycles = [None] * len(list(range(start, len(ts), fs)))
            for idx, frame in enumerate(range(start, len(ts), fs)):
                cycle_lengths = [len(tmp) for tmp in nx.cycle_basis(self.nxg[frame])]
                num_cycles[idx] = [cycle_lengths.count(i) for i in range(rl[0], rl[1])]
                a = str(ts[frame])
                b = " ".join([str(cycle_lengths.count(i)) for i in range(rl[0], rl[1])])
                c = str(len(cycle_lengths))
                tmp = fring.write(" ".join([a, b, c, "\n"]))

            # close file
            fring.close()
# End of class
def usage():
    """Usage:
            python3 ..\reaXtract.py -i bonds.reaxff.dump -r 3:10 -a 1:6,2:1,3:1,4:8,5:8,6:8,7:8,8:8
    """

if __name__ == "__main__":

    short_opts = "hi:b:a:r:"
    long_opts = ("--help",
                 "--in=","--infile=","--input=",
                 "--basename=","--base="
                 "--atom_map=","--atom_type_map=","--atm=",
                 "--rs=","--ring_size=","--ring_sizes=")

    try:
        optlist, args = getopt.getopt(sys.argv[1:],shortopts=short_opts,longopts=long_opts)
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)

    # Defaults:
    infile = ""
    basename = ""
    atom_type_map = ""
    eval_by = "type"
    count_rings = False
    ring_length = (3,10)

    # input options
    print(optlist)
    print(args)
    for o, a in optlist:
        if o == "-v":
            verbose = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-i", "--in=","--infile=","--input="):
            infile = a.strip()
        elif o in ("-b","--basename="):
            basename = a.strip()
        elif o in ("-a","--atom_map=","--atom_type_map=","--atm="):
            atom_type_map = a.strip()
        elif o in ("-e","--evaluate_by"):
            if a in ("type","element"):
                eval_by = a
            else: assert False,["unrecognized option: evaluate_by must be 'type' or 'element'",o,a]
        elif o in ("-r","--rings=","--countrings="):
            count_rings = True
            if len(a)>0:
                tmp = [int(i) for i in a.strip().split(":")]
                if len(tmp)==1: ring_length = (3,tmp[0])
                elif len(tmp)==2: ring_length = tuple(sorted([tmp[0],tmp[1]]))
            else:
                ring_length=ring_length
        else:
            assert False, ["unhandled option:",o,a]

    # check arguments
    assert os.path.isfile(infile),["input file not found:",infile]

    #rxt: ReaXtract = ReaXtract(infile=infile,
    #                           basename=basename,
    #                           atom_type_map="1:6,2:1,3:1,4:8,5:8,6:8,7:8,8:8",
    #                           ring_counter=count_rings, ring_limits=ring_length)
    rxt: ReaXtract = ReaXtract(basename=basename,
                               atom_type_map=atom_type_map)
    rxt.read(infile=infile)
    rxt.find_rxn()
    rxt.plot_rxn()
    if count_rings:
        print("ring_length",ring_length)
        rxt.count_rings(ring_limits=ring_length)
    print("finished")



