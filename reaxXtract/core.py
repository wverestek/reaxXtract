import os.path, sys, re
import random

from typing import TextIO, Union, List

import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from .utils import k_nearest_neighs, convert_str2dict
from .utils import ON2ELEM, ON2HEX, ELEM2HEX, DEFAULT_COLOR
from .reader import read_bonds
from .logger import log, configure_log

configure_log(level="DEBUG", force=True)

#import pdb
#pdb.set_trace()

class ReaxXtract:
    ##############
    # initialize #
    ##############
    def __init__(self, 
                 infile: Union[str, List[str]] = "", informat: str = "reaxff", basename: str = "",
                 atom_type_map: str = "",
                 checkframe: int = 1, stepframe: int = 1, stabiframes: int = 0,
                 hash_by: str = "type", 
                 rxn_bond_cutoff: int = 1, plot_bonds_cutoff: int = 5, seed: int = 42,
                 ring_counter: bool = False, ring_limits=(3, 10)):
        """
        A class to extract changes in bond topology over time.

        infile : str or list[str]
            A bond information file, a list of files, or glob pattern(s). The
            reader will accept a single filename or multiple files. For
            backward compatibility `self.infile` is the first filename (or "").
        infile : str or list[str]
            A bond information file, a list of files, or glob pattern(s). The
            reader will accept a single filename or multiple files. For
            backward compatibility `self.infile` is the first filename (or "").
        informat : str
            file type of the file containing bond information. Default "reaxff"
        basename : str
            base name for output. If not set the input file name is used as base name.
        startstep : int
            MD time step at start. Default: 0
        stopstep : int
            MD time step at end. Default: sys.maxsize (a very high umber)
        checkframe : int
            Number of frames difference to check for changed bonds. Default: 1
        stepframe : int
            Number of frames before the next evaluation is done. Default: 1
        hash_by : str
            hash for each reaction that is build upon 'element' or 'type'. 
            Hashes allow to identify if a similar reaction has already occured before or not.
        """
        
        log.info(f"Initializing ReaxXtract class object...")
        self.name: str = "ReaxXtract"

        # preserve the raw input(s) and provide a canonical list + first-file convenience attribute
        if isinstance(infile, (list, tuple)):
            self.infiles: List[str] = list(infile)
        elif infile is None or infile == "":
            self.infiles = []
        else:
            self.infiles = [infile]

        # backward-compatible single-file attribute (used in many places)
        self.infile: str = self.infiles[0] if len(self.infiles) > 0 else ""

        self.informat: str = informat.lower()
        if len(basename) > 0:
            self.basename:str = basename
        else:
            self.basename:str = re.sub(r'(\.(?:gz|txt|dat|data|dump))+$', '', os.path.basename(self.infile))
        
        #self.startstep:int = int(startstep)     # necessary?
        #self.stopstep:int = int(stopstep)       # necessary?
        #self.startidx:int = 0                   # necessary?
        #self.stopidx:int = 0                    # necessary?
        self.checkframe:int = int(checkframe)   # necessary?
        self.stabiframe:int = int(stabiframes)
        self.stepframe:int = int(stepframe)     # necessary?
        self.maxframe_offset = 0
        
        self.ts:list = []                                   # list of timesteps
        self.nxg:list = []                                  # list of networkx Graphs for each timestep
        self.frames:pd.DataFrame = pd.DataFrame(columns=("frame","timestep","graph")) # DataFrame with columns ['timestep','graph']

        self.rxns:pd.DataFrame = pd.DataFrame(
            columns=("frame","timestep","rxnID","rxnCount",
                     "edges_before","edges_after",
                     "atoms_rxn","atoms_env","atoms_plot",
                     "Gbefore","Gafter",
                     "rxn_hash_before","rxn_hash_after"))                          # DataFrame with reactions found
        #self.rxn_sets_before:list = [[]]                    # list of reaction sets before reaction
        #self.rxn_sets_after:list = [[]]
        #self.rxn_hashes_before:list = [[]]
        #self.rxn_hashes_after:list = [[]]
        self.hash_by:str = hash_by
        #self.rxn_elements_before = [[]]
        #self.rxn_elements_after = [[]]
        #self.rxn_id = []                                    # list of unique reaction IDs found
        #self.rxn_count = []                                 # count of each reaction found

        self.ring_counter = ring_counter
        self.ring_limits = ring_limits
        self.ring_sets = []

        self.rxn_bond_cutoff:int = int(rxn_bond_cutoff)     # number of bonds distance that are considered one reaction, default: 0
        self.plot_bonds_cutoff:int = int(plot_bonds_cutoff)               # number of bonds distance to include for plots
        random.seed(a=int(seed))                            # reproducibility of pyplot plots
        
        # expose ON2ELEM/ON2HEX as instance attributes so the rest of the class can use them
        self.on2elem:dict = ON2ELEM
        self.on2hex:dict = ON2HEX
        self.elem2hex:dict = ELEM2HEX
        self.default_color:str = DEFAULT_COLOR
        # atom mapping from type to ordinal number -> user input as str
        self.atom_type_map = atom_type_map
        self.type2on:dict = convert_str2dict(self.atom_type_map)
        #self.elem2hex:dict = {self.on2elem[i]:self.on2hex[i] for i in range(1,len(self.on2elem)+1)}
        self.elem2hex:dict = {self.on2elem[i]:self.on2hex[i] for i in self.on2elem.keys()}

    ##############
    # read bonds #
    ##############
    def read(self, infile:str=None, informat:str=None):
        """
        Read bond file and populate instance state:
          - self.ts, self.nxg
          - self.startidx / self.stopidx (based on startstep/stopstep)
          - self.frames (pd.DataFrame) with columns ['timestep','graph']
        Returns the frames DataFrame.
        """
        # choose infile to pass to reader.read_bonds:
        # - if caller provided `infile` use that
        # - else prefer self.infiles (list) if present, otherwise self.infile (string)
        if infile is not None:
            infile_local = infile
        else:
            infile_local = self.infiles if len(self.infiles) > 0 else self.infile

        informat_local = informat or self.informat
        if not infile_local:
            raise ValueError("No infile supplied to read()")

        # read bond file(s)
        ts, nxg = read_bonds(infile_local, informat_local)
        # set element attribute for each node in each graph
        for i,g in enumerate(nxg):
            for n, data in g.nodes(data=True):
                atom_type = data.get("type", None)
                if atom_type is not None and atom_type in self.type2on:
                    element = self.type2on[atom_type]
                    nx.set_node_attributes(nxg[i], {n: element}, name="element")
                else:
                    nx.set_node_attributes(nxg[i], {n: element}, name="element")
                    log.warn(f"Warning: atom type {atom_type} not in atom_type_mapping, setting element to 'X'",stacklevel=2)
        
        frames_arr = list(range(len(ts)))
        if self.frames.empty:
            # first read, no existing frames, just set frames DataFrame
            self.frames = pd.DataFrame({"frame":frames_arr,"timestep": ts,"graph": nxg})
        else:
            # subsequent read
            # drop frames in case large files are read in multiple calls to avoid 
            # memory issues, but keep stabilize frames for reaction checking
            maxframe_keep =  self.frames["frame"].iloc[-1] - self.stabiframe
            idx = np.where(self.frames["frame"].lt(maxframe_keep))[0].tolist()
            _ = self.frames.drop(index=idx, inplace=True)
            # renumber new frames to continue from last frame + 1
            maxframe_old = self.frames["frame"].max() if not self.frames.empty else -1
            frames_arr = [i + maxframe_old + 1 for i in frames_arr]
            # concatenate new frames to existing frames DataFrame
            self.frames = pd.concat([self.frames, 
                                     pd.DataFrame({"frame":frames_arr,"timestep": ts,"graph": nxg})],
                                    ignore_index=True)

        # determine start and stop indices
        ### necessary?
        #if self.startstep in self.frames:
        #    self.startidx = ts.index(self.startstep)
        #else:
        #    self.startidx = 0
        #    log.warn(f"Warning: startstep {self.startstep} not found in timesteps, setting to first timestep {ts[self.startidx]}")    
        #if self.stopidx > len(ts) - 1 or self.stopstep not in ts:
        #    self.stopidx = len(ts) - 1
        #    log.warn(f"Warning: stopstep {self.stopstep} not found in timesteps, setting to last timestep {ts[self.stopidx]}")
        #else:
        #    self.stopidx = ts.index(self.stopstep)




    ##################
    # find reactions #
    ##################
    def find_rxns(self):
        """
        read_rxns(self,
                  infile: Optional[str] = None,
                  informat: Optional[str] = None,
                  force: bool = False) -> pandas.DataFrame

        Convenience method that runs reaction detection, 
        returning a reaction-level pandas.DataFrame.

        Summary
        - Calls `self.read(infile, informat)` (unless already read and force is False)
          and then `self.find_rxns()`.
        - After execution the instance contains:
            - self.frames (DataFrame, one row per frame)
            - self.rxns (DataFrame, one row per detected reaction occurrence)
            - self.rxn_id and self.rxn_count for unique-reaction bookkeeping

        Returns
        - pandas.DataFrame: `self.rxns` - reaction occurrences with typical columns:
            `frame_idx`, `timestep`, `rxn_idx`, `rxn_id`, `rxn_hash`,  `rxn_total_count`,
            `before_components`, `after_components`, `types_before`, `types_after`,
            `elements_before`, `elements_after`, `atoms_list`.

        Side effects
        - Overwrites/sets: `self.rxns`, `self.rxn_id`, `self.rxn_count`.
        - Writes summary file `self.basename + "_rxnIDs.dat"` (same behavior as current find_rxns).

        Errors / Exceptions
        - ValueError: when no input file available or start/stop steps invalid.
        - FileNotFoundError or reader errors propagated from `read_bonds`.
        - RuntimeError for unexpected internal state.

        Performance notes
        - WL hashing and graph operations can be CPU / memory intensive for large trajectories.
          Consider chunked processing, serializing graphs to disk, or disabling hashing for
          very large datasets. For chunked processing rxn_ids might be inconsistent across chunks.

        Example
            r = ReaxXtract(infile="bonds.reaxff.dump")
            rxns = r.read_rxns()           # reads frames and finds reactions
            rxns.to_csv("rxns_summary.csv")
        """

        if self.frames.empty:
            log.info(f"No bond data read yet, reading now... {self.infile}")
            self.read()

        log.info("Searching reactions...")
        f_rxn:TextIO = open(self.basename + "_rxnIDs.dat", "wt")
        mystr = "# Timestep \t RxnID \t RxnCount \t FromIDs:ToIDs \t FromType:ToType \t FromElem:ToElem \t Rxn_hashes"
        print(mystr)
        f_rxn.write(mystr+"\n")
        

        start = 0 
        stop =  self.frames["frame"].count() - self.stabiframe 
        cf = self.checkframe
        fs = self.stepframe

        #self.rxn_id = []
        #self.rxn_count = []

        df_file = pd.DataFrame(columns=self.rxns.columns)
        
        for idx in range(start + cf, stop, fs):
            # dicts for conversion
            node2element = nx.get_node_attributes(self.frames["graph"].iloc[idx], name="element")
            node2type    = nx.get_node_attributes(self.frames["graph"].iloc[idx], name="type")

            before_idx = idx - cf
            after_idx = idx

            Gbefore = self.frames["graph"].iloc[idx-cf]
            Gafter = self.frames["graph"].iloc[idx]

            log.info(f"evaluating frame: {idx} ({self.frames["frame"].iloc[idx]}) timestep: {self.frames['timestep'].iloc[idx]}")
            # find reacting atoms
            [edges_sets_before,edges_sets_after, reaction_sets] = self.find_reacting_atoms_for_two_frames(Gbefore, Gafter)   # list(set(int,),set(int,))
            log.debug(f"reaction_sets: {reaction_sets} with edges_before: {edges_sets_before} and edges_after: {edges_sets_after}")
                
            if len(reaction_sets) > 0:
                df_frame = self.rsets_to_pd(before_idx, after_idx,
                                            reaction_sets, edges_sets_before, edges_sets_after)
                if df_frame is not None and not df_frame.empty:
                    df_file = pd.concat([df_file, df_frame], ignore_index=True)
                else:
                    log.warn("You should not be here. Maybe you have discovered a bug. Please consider reporting with a minimal example")

        f_rxn.close()
        
        # adding newly found reactions to self.rxns DataFrame
        # first search
        if self.rxns.empty:
            self.rxns = df_file.copy()
        else:
            log.info(f"appending new search to already stored reactions")
            df_file["frame"] = df_file["frame"] + self.maxframe_offset
            self.rxns = pd.concat([self.rxns,df_file],ignore_index=True)

        #self.df1 = df_file.copy()
        # renumber reactions and count unique reactions
        self.renumber_and_count_rxns()


    # find reacting atoms for two frames #
    def find_reacting_atoms_for_two_frames(self,Gbefore:nx.Graph,Gafter:nx.Graph):
        # compare graphs, edges that are not in both graphs => reaction
        #reacting_edges = nx.symmetric_difference(Gbefore,Gafter).edges()
        reacting_edges_before = nx.difference(Gbefore,Gafter).edges()
        reacting_edges_after  = nx.difference(Gafter,Gbefore).edges()
        reacting_edges = set(reacting_edges_before).union(set(reacting_edges_after))
        log.debug(f"reacting_edges: {reacting_edges}")
        # atom IDs that are involved with changed bond connectivity
        # reacting_atoms = set(i for j in reacting_edges for i in j)          # set(int)
        reacting_atoms = set(i for j in set(reacting_edges) for i in j)          # set(int)
        log.debug(f"reacting_atoms: {reacting_atoms}")
        
        
        # reactions found
        if len(reacting_atoms) > 0:
            # include atoms rxn_bond_cutoff bonds away
            if self.rxn_bond_cutoff > 0:
                reacting_atoms1 = k_nearest_neighs(Gbefore,reacting_atoms,self.rxn_bond_cutoff)
                reacting_atoms2 = k_nearest_neighs(Gafter ,reacting_atoms,self.rxn_bond_cutoff)  
                reacting_atoms = reacting_atoms1.union(reacting_atoms2)     # combine sets
            # group reacting atoms by connectivity to find individual        
            # create small graph with reacting atoms and edges (before and after)
            Gcombine = Gbefore.subgraph(reacting_atoms).copy()              # nx.Graph
            Gcombine.add_edges_from(reacting_edges)                         # add reacting edges 
            # group combined edges (before and after) by connectivity 
            # reactions are distinguished by connected groups of edges
            # list of sets, one set per found reaction
            reacting_atoms_sets = list(nx.connected_components(Gcombine))   # list(set(int,),)
            log.debug(f"reacting_atoms_sets: {reacting_atoms_sets}")

            # check for overlapping reaction cores in reaction sets
            reacting_edges_before_sets = list(set())
            reacting_edges_after_sets = list(set())
            for i,rset in enumerate(reacting_atoms_sets):
                reacting_edges_before_sets.append([set(bond) for bond in reacting_edges_before if set(bond).issubset(rset)])
                reacting_edges_after_sets.append( [set(bond) for bond in reacting_edges_after  if set(bond).issubset(rset)])
            log.debug(f"reacting_edges_before {reacting_edges_before}")
            log.debug(f"reacting_edges_after {reacting_edges_after}")
            
        else:
            reacting_edges_before_sets = []
            reacting_edges_after_sets = []
            reacting_atoms_sets = []

        return reacting_edges_before_sets, reacting_edges_after_sets, reacting_atoms_sets


    # convert reaction sets to pandas DataFrame format for further analysis and plotting #
    def rsets_to_pd(self,before:int,after:int,reaction_sets,edges_sets_before,edges_sets_after):
        """
        Convert reaction sets and their associated edge changes into 
        a pandas DataFrame format for further analysis and plotting.
        Parameters:
        before (int): Index of the "before" frame.
        after (int): Index of the "after" frame.
        reaction_sets (list of sets): List of sets of atom IDs involved in each reaction.
        edges_sets_before (list of sets): edges that vanish in this reaction.
        edges_sets_after (list of sets): edges that become created in this reaction.
        """
        tmp_list = []
        # before and after sets for each individual reaction
        for i,(rset,edges_before,edges_after) in enumerate(zip(reaction_sets,edges_sets_before,edges_sets_after)):
            log.debug(f"rset: {rset} with edges_before {edges_before} and edges_after {edges_after}")
            log.debug(f"rset: edge atoms {list(set().union(*edges_before).union(*edges_after))}")
            atoms_rxn = list(set().union(*edges_before).union(*edges_after))
            atoms_env = list(rset)
            plot_atoms1 = k_nearest_neighs(self.frames["graph"].iloc[before], atoms_env, self.plot_bonds_cutoff)
            plot_atoms2 = k_nearest_neighs(self.frames["graph"].iloc[after] , atoms_env, self.plot_bonds_cutoff)  
            atoms_plot = list(plot_atoms1.union(plot_atoms2))     # combine sets
            
            # before including plot atoms 
            Gbefore = self.frames["graph"].iloc[before].subgraph(atoms_plot).copy()
            hash_before = nx.weisfeiler_lehman_graph_hash(Gbefore.subgraph(atoms_env), node_attr=self.hash_by)
            # after including plot atoms 
            Gafter = self.frames["graph"].iloc[after].subgraph(atoms_plot).copy()
            hash_after = nx.weisfeiler_lehman_graph_hash(Gafter.subgraph(atoms_env), node_attr=self.hash_by)

            tmp_list.append({"frame":self.frames["frame"].iloc[after],
                             "timestep":self.frames["timestep"].iloc[after],
                             "rxnID":None,"rxnCount":None,
                             "edges_before":edges_before,"edges_after":edges_after,
                             "atoms_rxn":sorted(atoms_rxn),"atoms_env":sorted(atoms_env),"atoms_plot":sorted(atoms_plot),
                             "Gbefore":Gbefore,"Gafter":Gafter,
                             "rxn_hash_before":hash_before,"rxn_hash_after":hash_after})
        return pd.DataFrame(tmp_list)

    
    # renumber reactions and count unique reactions #
    def renumber_and_count_rxns(self):
        if self.rxns.empty:
            log.info("No reactions to renumber and count.")
            return

        crxns = self.rxns["rxn_hash_before"].count() # total count of reactions
        rxn_hashes = (self.rxns["rxn_hash_before"] + [":"]*crxns + self.rxns["rxn_hash_after"]).tolist()
        hash2id = {h: i for i, h in enumerate(set(rxn_hashes))}     # dict: rxn_hash => rxn_id
        nrxns = np.max(list(hash2id.values()))+1                    # number of individual reactions

        rxn_id = np.array([None] * crxns)
        rxn_count = np.array([0] * crxns)
        counter_arr = np.array([0] * nrxns)

        for idx, h in enumerate(rxn_hashes):
            rxn_id[idx] = hash2id[h]
            counter_arr[hash2id[h]] += 1
            rxn_count[idx] = counter_arr[hash2id[h]]
        
        self.rxns["rxnID"] = rxn_id
        self.rxns["rxnCount"] = rxn_count
        

    # remove reactions that reverse in stabilize frames #
    def remove_reversing_reactions(self,nframes:int=None,df:pd.core.frame.DataFrame=None):
        self.stabiframe = nframes or self.stabiframe
        
        # use on self.rxns or another Pandas DataFrame if provided as argument
        df = df or self.rxns

        rmv_idx = []
        for idx,row in self.rxns.iterrows():
            hash_before = row["rxn_hash_before"]
            hash_after = row["rxn_hash_after"]
            current_frame = row["frame"]
            max_frame = current_frame + self.stabiframe
            
            mask_frame = df["frame"].gt(current_frame) & df["frame"].le(max_frame)
            mask_hash = (df["rxn_hash_before"] == hash_after) & (df["rxn_hash_after"] == hash_before)
            # list comparison: convert to tuple
            target_atoms = tuple(row["atoms_env"])
            mask_atoms = df["atoms_env"].map(tuple) == target_atoms
            mask = mask_frame & mask_hash & mask_atoms
            tmp = np.where(mask)[0]
            if len(tmp) == 1:
                rev_idx = tmp[0]
                rmv_idx.append(idx)
                rmv_idx.append(rev_idx)
        
        if len(rmv_idx) > 0:
            log.info(f"Reaction(s) found that reverses within {self.stabiframe} frames, removing reactions")
            #log.info(f"{df.iloc[rmv_idx]}")
            #df_rmv = df.iloc[rmv_idx].copy()
            df_rmv = df.drop(rmv_idx, inplace=True)
            df.reset_index(drop=True, inplace=True)
            self.renumber_and_count_rxns()
        else:
            df_rmv = None
        
        return df_rmv

    # plot reactions #
    def plot_rxns(self):
        if self.rxns.empty:
            log.warn("No reactions found to plot.")
            return
        cf = self.checkframe
        outfolder = self.basename + ".dir_df"
        if not os.path.exists(outfolder):
            os.makedirs(outfolder, exist_ok=True)

        for idx, rxn in self.rxns.iterrows():
            frame = rxn["frame"]
            timestep = rxn["timestep"]
            
            Gbefore = rxn["Gbefore"]
            Gafter = rxn["Gafter"]
            #active_edges_before = nx.difference(Gbefore, Gafter).edges()
            #active_edges_after  = nx.difference(Gafter, Gbefore).edges()
            active_edges = nx.symmetric_difference(Gbefore, Gafter).edges()

            # initial positions
            Gprint = Gbefore.copy()
            Gprint.add_edges_from(active_edges)
            pos = nx.kamada_kawai_layout(Gprint)
            pos = nx.spring_layout(Gprint, iterations=200, pos=pos)
            del Gprint

            
            plt.clf()
            # before reaction
            # prepare data
            plt.figure(figsize=(14,6))
            plt.suptitle(f"Frame: {frame} Timestep: {timestep}")
            plt.subplot(1, 2, 1)
            plt.title("Before")
            node2elem = nx.get_node_attributes(Gbefore, name="element")
            node2type = nx.get_node_attributes(Gbefore, name="type")
            color = [self.elem2hex.get(v, self.default_color) for v in
                      nx.get_node_attributes(Gbefore, name="element").values()]
            node_labels = {k:  str(node2elem[k]) + ":" + str(node2type[k]) + "\n" + str(k) for k in Gbefore}
            if len(nx.get_edge_attributes(Gbefore,"bo")) > 0:
                bo = [tmp for tmp in list(nx.get_edge_attributes(Gbefore, "bo").values())]
            else:
                bo = [1.0 for i in Gbefore.edges()]
            pos_before = nx.spring_layout(Gbefore, iterations=50, pos=pos)
            # plot data
            #nx.draw_networkx_edges(Gbefore, edgelist=active_edges_before, alpha=0.6, width=5.0, edge_color="tab:red", pos=pos)
            nx.draw_networkx_edges(Gbefore, edgelist=active_edges, alpha=0.6, width=5.0, edge_color="tab:red", pos=pos_before)
            nx.draw(Gbefore, pos=pos_before, node_color=color, 
                    with_labels=True, labels=node_labels, font_size=6,
                    node_size=300, edge_color="black", width=bo)
            plt.tight_layout()
            
            # after reaction
            # prepare data
            plt.subplot(1, 2, 2)
            plt.title("After")
            color = [self.elem2hex.get(v, self.default_color) for v in
                     nx.get_node_attributes(Gafter, name="element").values()]
            node_labels = {k:  str(node2elem[k]) + ":" + str(node2type[k]) + "\n" + str(k) for k in Gafter}
            if len(nx.get_edge_attributes(Gafter,"bo")) > 0:
                bo = [tmp for tmp in list(nx.get_edge_attributes(Gbefore, "bo").values())]
            else:
                bo = [1.0 for i in Gbefore.edges()]
            pos_after = nx.spring_layout(Gafter, iterations=50, pos=pos)
            # plot data
            #nx.draw_networkx_edges(Gafter, edgelist=active_edges_after, alpha=0.6, width=5.0, edge_color="tab:red", pos=pos)
            nx.draw_networkx_edges(Gafter, edgelist=active_edges, alpha=0.6, width=5.0, edge_color="tab:red", pos=pos_after)
            nx.draw(Gafter, pos=pos_after, node_color=color, 
                    with_labels=True, labels=node_labels, font_size=6,
                    node_size=300, edge_color="black", width=bo)
            plt.tight_layout()
            f_out = os.path.join(outfolder, self.basename + "_" + str(timestep) + "_Rxn" + str(idx) + ".png")
            plt.savefig(f_out,dpi=200)


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
