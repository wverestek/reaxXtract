#from fileinput import filename
from math import e
import os, os.path

from networkx.algorithms.operators import union
from networkx.drawing import draw
os.environ.setdefault("MPLBACKEND", "Agg")
import sys, re
import random

from typing import TextIO, Union, List

import matplotlib
# Use non-interactive backend to avoid Qt/X11 errors in headless environments
matplotlib.use("Agg")

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
                 ring_counter: bool = False, loop_limits:tuple[int,int]=None):
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

        # store input(s): `infiles` may be a single filename (str) or a list of filenames
        # reader.read_bonds accepts either a string, a list or glob pattern(s).
        self.infile: Union[str, List[str]] = infile
        self.informat: str = informat.lower()
        # derive basename based on only or first filename string
        if isinstance(self.infile, str) and len(self.infile) > 0:
            self.basename: str = re.sub(r'(\.(?:gz|txt|dat|data|dump))+$', '', os.path.basename(self.infile))
        elif isinstance(self.infile, list):
            self.basename: str = re.sub(r'(\.(?:gz|txt|dat|data|dump))+$', '', os.path.basename(self.infile[0]))
        
        self.checkframe:int = int(checkframe)   # necessary?
        self.stabiframe:int = int(stabiframes)
        self.stepframe:int = int(stepframe)     # necessary?
        self.maxframe_offset = 0
        
        # pandas DataFrame to store global bond topology
        self.frames:pd.DataFrame = pd.DataFrame(columns=("frame","timestep","graph")) # DataFrame with columns ['timestep','graph']
        
        # pandas DataFrame to store reactions found
        self.rxns:pd.DataFrame = pd.DataFrame(
            columns=("frame","timestep","rxnID","rxnCount",
                     "edges_before","edges_after",
                     "atoms_rxn","atoms_env","atoms_plot",
                     "Gbefore","Gafter",
                     "rxn_hash_before","rxn_hash_after"))                          # DataFrame with reactions found
        
        # build hash based on node attribute 'element' or 'type'
        self.hash_by:str = hash_by                      
        
        self.rxn_bond_cutoff:int = int(rxn_bond_cutoff)         # number of bonds distance that are considered one reaction, default: 0
        self.plot_bonds_cutoff:int = int(plot_bonds_cutoff)     # number of bonds distance to include for plots
        random.seed(a=int(seed))                                # reproducibility of pyplot plots

        # find and count loops in each global frame        
        self.loop_limits = loop_limits

        # expose ON2ELEM/ON2HEX as instance attributes so the rest of the class can use them
        self.on2elem:dict = ON2ELEM
        self.on2hex:dict = ON2HEX
        self.elem2hex:dict = ELEM2HEX
        self.default_color:str = DEFAULT_COLOR
        # atom mapping from type to ordinal number -> user input as str
        self.atom_type_map = atom_type_map
        self.type2on:dict = convert_str2dict(self.atom_type_map)
        self.elem2hex:dict = {self.on2elem[i]:self.on2hex[i] for i in self.on2elem.keys()}

    # read bonds #
    def read(self, infile: Union[str, List[str]] = None, informat: str = None):
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
            # use self.infile (string or list)
            infile_local = self.infile

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
                    nx.set_node_attributes(nxg[i], {n: "X"}, name="element")
                    log.warn(f"Warning: atom type {atom_type} not in atom_type_mapping, setting element to 'X'",stacklevel=2)
        
        frames_arr = list(range(len(ts)))
        if self.frames.empty:
            # first read, no existing frames, just set frames DataFrame
            self.frames = pd.DataFrame({"frame":frames_arr,"timestep": ts,"graph": nxg})
        else:
            # subsequent read
            # drop frames in case large files are read in multiple calls to avoid 
            # memory issues, but keep stabilize frames for reaction checking
            nframes = len(self.frames["frame"])
            if  nframes > 100 and nframes > self.stabiframe:
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
        # check timesteps are monotonic increasing
        if not self.frames["timestep"].is_monotonic_increasing:
            log.warning(f"READER: Timesteps are not monotonic increasing!!!")


    # find reactions #
    def find_reactions(self):
        """
        read_rxns(self,
                  infile: Optional[str] = None,
                  informat: Optional[str] = None,
                  force: bool = False) -> pandas.DataFrame

        Convenience method that runs reaction detection, 
        returning a reaction-level pandas.DataFrame.

        Summary
        - Calls `self.read(infile, informat)` (unless already read and force is False)
          and then `self.find_reactions()`.
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
        
        # check for ascending timesteps
        if not self.frames["timestep"].is_monotonic_increasing:
            log.warning(f"Timesteps are not monotonic increasing!!!")

        log.info("Searching reactions...")
        f_rxn:TextIO = open(self.basename + "_rxnIDs.dat", "wt")
        mystr = "# Timestep \t RxnID \t RxnCount \t FromIDs:ToIDs \t FromType:ToType \t FromElem:ToElem \t Rxn_hashes"
        print(mystr)
        f_rxn.write(mystr+"\n")
        

        start = 0 
        stop =  len(self.frames["frame"]) - self.stabiframe 
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

            log.info(f"evaluating frame: {idx} ({self.frames['frame'].iloc[idx]}) timestep: {self.frames['timestep'].iloc[idx]}")
            # find reacting atoms
            [edges_sets_before,edges_sets_after, reaction_sets] = self._find_reacting_atoms_for_two_frames(Gbefore, Gafter)   # list(set(int,),set(int,))
            log.debug(f"reaction_sets: {reaction_sets} with edges_before: {edges_sets_before} and edges_after: {edges_sets_after}")
                
            if len(reaction_sets) > 0:
                df_frame = self._rsets_to_pd(before_idx, after_idx,
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
        renumber_and_count_rxns()


    # find reacting atoms for two frames #
    def _find_reacting_atoms_for_two_frames(self,Gbefore:nx.Graph,Gafter:nx.Graph):
        # compare graphs, edges that are not in both graphs => reaction
        #reacting_edges = nx.symmetric_difference(Gbefore,Gafter).edges()
        all_edges_before = nx.difference(Gbefore,Gafter).edges()
        all_edges_after  = nx.difference(Gafter,Gbefore).edges()
        all_reacting_edges = set(all_edges_before).union(set(all_edges_after))
        log.debug(f"all_reacting_edges: {all_reacting_edges}")
        # atom IDs that are involved with changed bond connectivity
        reacting_atoms_core = set(i for j in set(all_reacting_edges) for i in j)       # set(int,)
        log.debug(f"reacting_atoms_core: {reacting_atoms_core}")

        # subgraph with reacting atoms and edges (before and after) 
        # to find connected components as individual reactions
        Gcombine = Gbefore.subgraph(reacting_atoms_core).copy()              # nx.Graph
        Gcombine.add_edges_from(all_reacting_edges)                         # nx.Graph, add reacting edges 
        reacting_atoms_sets = list(nx.connected_components(Gcombine))       # list(set(int,),)
        
        nsets = len(reacting_atoms_sets)
        # reactions found
        if nsets > 0:
            
            # include atoms rxn_bond_cutoff bonds away
            if self.rxn_bond_cutoff > 0:
                
                # expand individual sets by rxn_bond_cutoff bonds
                log.debug(f"reacting_atoms_sets before expansion: {reacting_atoms_sets}")
                tmpsets = []
                for i,iset in enumerate(reacting_atoms_sets):
                    reacting_atoms1 = k_nearest_neighs(Gbefore,iset,self.rxn_bond_cutoff)
                    reacting_atoms2 = k_nearest_neighs(Gafter ,iset,self.rxn_bond_cutoff)  
                    #tmpsets[i] = reacting_atoms1.union(reacting_atoms2)     # combine sets
                    tmpsets.append(reacting_atoms1.union(reacting_atoms2))    # combine sets
                log.debug(f"reacting_atoms_sets after expansion: {tmpsets}")

                # more than one set, check if sets have common atoms after expansion
                # and merge if necessary, otherwise just use the expanded sets as reaction sets
                if len(tmpsets) > 1:
                    log.debug(f"Merging:")
                    reacting_atoms_sets = []
                    # merge sets that have common atoms after expansion
                    log.debug(f"reacting_atoms_sets before merging: {reacting_atoms_sets}")
                    for i,iset in enumerate(tmpsets):
                        for j in range(i+1, nsets):
                            if len(iset.intersection(tmpsets[j])) > 0:
                                reacting_atoms_sets.append(iset.union(tmpsets[j]))
                                tmpsets[j] = set()    # reomve from further consideration
                    log.debug(f"reacting_atoms_sets after merging: {reacting_atoms_sets}")

                else:
                    # only one set, just use the expanded set
                    reacting_atoms_sets = tmpsets
                    
            # no expansion, just use original sets of reacting atoms and edges
            else:
                pass
            log.debug(f"reacting_atoms_sets: {reacting_atoms_sets}")

            # construct edge sets for each reaction set
            reacting_edges_before_sets = list(set())
            reacting_edges_after_sets = list(set())
            for i,rset in enumerate(reacting_atoms_sets):
                reacting_edges_before_sets.append([set(bond) for bond in all_edges_before if set(bond).issubset(rset)])
                reacting_edges_after_sets.append( [set(bond) for bond in all_edges_after  if set(bond).issubset(rset)])
            log.debug(f"reacting_edges_before: {reacting_edges_before_sets}")
            log.debug(f"reacting_edges_after: {reacting_edges_after_sets}")
        
        # no reactions found, return empty lists
        else:
            reacting_edges_before_sets = []
            reacting_edges_after_sets = []
            reacting_atoms_sets = []

        return reacting_edges_before_sets, reacting_edges_after_sets, reacting_atoms_sets


    # convert reaction sets to pandas DataFrame format for further analysis and plotting #
    def _rsets_to_pd(self,before:int,after:int,reaction_sets,edges_sets_before,edges_sets_after):
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
# End of class

## work on reactions and topology ##
# renumber reactions and count unique reactions #
def renumber_and_count_rxns(self, df:pd.core.frame.DataFrame=None) -> pd.core.frame.DataFrame:
    """
    Renumber reactions and count unique reactions based on their hashes in pandas DataFrame. 
    This method updates two columns of the DataFrame: 'rxnID' and 'rxnCount'.
    If df is None, operates on self.rxns and updates it in place. Otherwise, operates on the 
    provided DataFrame and returns a new DataFrame with the updated columns.
    """
    # choose target DataFrame
    df_in = self.rxns if df is None else df
    # operate on a copy to avoid surprising in-place side effects for caller
    df_work = df_in.copy()

    if df_work is None or df_work.empty:
        log.info("No reactions to renumber and count.")
        return

    # reset index to ensure consistent indexing for reaction ID assignment
    df_work.reset_index(drop=True, inplace=True)

    crxns = len(df_work["rxn_hash_before"]) # total count of reactions
    rxn_hashes = df_work["rxn_hash_before"] + [":"]*crxns + df_work["rxn_hash_after"]
    id2hash = dict(enumerate(pd.unique(rxn_hashes))).items()    # unique rxn_hashes in order of appearance
    hash2id = dict((v,k) for k,v in id2hash)                    # unique rxn_hashes in order of appearance
    nrxns = len(hash2id.keys())                                 # number of individual reactions

    rxn_id = np.array([None] * crxns)
    rxn_count = np.array([0] * crxns)
    counter_arr = np.array([0] * nrxns)

    for idx, h in enumerate(rxn_hashes):
        rxn_id[idx] = hash2id[h]
        counter_arr[hash2id[h]] += 1
        rxn_count[idx] = counter_arr[hash2id[h]]
    
    df_work["rxnID"] = rxn_id
    df_work["rxnCount"] = rxn_count

    return df_work


# remove reactions that reverse in stabilize frames #
def filter_transient_reactions(self, nframes:int=None, df:pd.core.frame.DataFrame=None)-> list[pd.core.frame.DataFrame,pd.core.frame.DataFrame]:
    """
    Remove reactions that reverse within a certain number of frames (stabiframe) after their 
    occurrence. Based on 'frame', 'rxn_hash_before' and 'rxn_hash_after' in pandas DataFrame.
    Returns a tuple of two DataFrames: (filtered_reactions, removed_reactions) where:
    - filtered_reactions: DataFrame with reactions that do not reverse within stabiframe frames.
    - removed_reactions: DataFrame with reactions that reverse within stabiframe frames.
    """
    # choose target DataFrame
    df_in = self.rxns if df is None else df
    # operate on a copy to avoid surprising in-place side effects for caller
    df_work = df_in.copy()

    self.stabiframe = nframes or self.stabiframe

    rmv_idx = []
    for idx,row in df_work.iterrows():
        hash_before = row["rxn_hash_before"]
        hash_after = row["rxn_hash_after"]
        current_frame = row["frame"]
        max_frame = current_frame + self.stabiframe
        
        mask_frame = df_work["frame"].gt(current_frame) & df_work["frame"].le(max_frame)
        mask_hash = (df_work["rxn_hash_before"] == hash_after) & (df_work["rxn_hash_after"] == hash_before)
        # list comparison: convert to tuple
        target_atoms = tuple(row["atoms_env"])
        mask_atoms = df_work["atoms_env"].map(tuple) == target_atoms
        mask = mask_frame & mask_hash & mask_atoms
        tmp = np.where(mask)[0]
        if len(tmp) > 0:
            rev_idx = tmp[0]
            rmv_idx.append(idx)
            rmv_idx.append(rev_idx)
    
    if len(rmv_idx) > 0:
        log.info(f"{len(rmv_idx)} Reaction(s) found that reverse within {self.stabiframe} frames, removing reactions")
        df_rmv = df_work.drop(rmv_idx, inplace=True)
        df_work = renumber_and_count_rxns(df_work)
        df_rmv  = renumber_and_count_rxns(df_rmv)
        renumber_and_count_rxns()
    else:
        df_work = renumber_and_count_rxns(df_work)
        df_rmv = None

    return df_work, df_rmv


# remove_atoms #
def remove_atoms_by_type(df:pd.core.frame.DataFrame=None, target_atoms:tuple[int|str,...]=None) -> pd.core.frame.DataFrame:
    # copy DataFrame to avoid surprising in-place side effects for caller
    df_work = df.copy(deep=True)
            
    for idx, row in df_work.iterrows():
        # real copy to avoid modifying the original graph in self.frames
        nxg = row["graph"]

        if target_atoms is None:
            nodes = list(nxg.nodes())
        else:
            nodes = [node for node, node_data in nxg.nodes(data=True)
                     if node_data.get('type') in target_atoms 
                     or node_data.get('element') in target_atoms]

        nxg.remove_nodes_from(nodes)
        
        # update the graph in the DataFrame with the modified graph
        # should be unnecessary since we are modifying the graph in place, but to be explicit:
        df_work.at[idx, "graph"] = nxg
    
    # return independent DataFrame with modified graphs
    return df_work


# remove atoms by somorph search of subgraph #
def remove_atoms_by_pattern(template_node_ids:list|set, delete_node_ids:list|set, df:pd.core.frame.DataFrame=None, node_attr:str='type', pattern_from_frame:int=0) -> pd.core.frame.DataFrame:
    """
    Simplifies the graph by matching a template pattern and removing specific 
    nodes. Handles molecular symmetry by filtering unique node sets.
    """
    from networkx.algorithms import isomorphism

    # 1. Sanity Check
    template_set = set(template_node_ids)
    delete_set = set(delete_node_ids)
    if not delete_set.issubset(template_set):
        invalid_ids = list(delete_set - template_set)
        log.error(f"Nodes {invalid_ids} are not in template_node_ids!")
        raise ValueError("Invalid delete_node_ids provided.")
    
    df_work = df.copy(deep=True)

    # 2. Create the template graph
    template = df_work["graph"].iloc[pattern_from_frame].subgraph(template_node_ids).copy()
    log.info(f"Starting topology reduction")
    log.info(f"Template pattern nodes: {list(template.nodes())}")

    nm = nx.isomorphism.categorical_node_match(node_attr, None)

    for idx, (df_idx, frame) in enumerate(df_work.iterrows()):
        nxg = frame["graph"]
        ncomp_before = nx.number_connected_components(nxg)

        # 3. Setup the GraphMatcher
        gm = nx.isomorphism.GraphMatcher(nxg, template, node_match=nm)

        # 4. Identify UNIQUE matches (Filtering Symmetries)
        nodes_to_remove = set()
        unique_molecule_footprints = set()

        for match in gm.subgraph_isomorphisms_iter():
            # match.keys() found pattern ID, match.values() corresponding template ID
            match_dict = dict(match)
            molecule_footprint = frozenset(match_dict.keys())
        
            if molecule_footprint not in unique_molecule_footprints:
                unique_molecule_footprints.add(molecule_footprint)
                # add key (pattern id) if value (template id) is in delete_node_ids
                #nodes_to_remove.add([k for k,v in match_dict.items() if v in delete_node_ids])
                ## Nur für den ersten gefundenen Isomorphismus dieses Moleküls löschen
                inv_match = {v: k for k, v in match.items()}
                for d_id in delete_node_ids:
                    if d_id in inv_match:
                        nodes_to_remove.add(inv_match[d_id])
    
        match_count = len(unique_molecule_footprints)

        if match_count == 0:
            log.warning(f"Index {idx}: No matches found! Check atom types.")
            continue

        # 5. Global removal
        initial_count = nxg.number_of_nodes()
        nxg.remove_nodes_from(nodes_to_remove)
    
        # 6. Connectivity Check
        ncomp_after = nx.number_connected_components(nxg)
        log.info(f"Index {idx}: Pattern found {match_count} times. Reduced nodes from {initial_count} to {nxg.number_of_nodes()}.")

        if ncomp_after > ncomp_before:
            log.warning(f"Index {idx}: Network is no longer fully connected! Fragments: {ncomp_before} vs. {ncomp_after}")
    
    return df_work


# plot reactions #
def plot_rxns(self, df:pd.core.frame.DataFrame=None, basename:str=None) -> None:
        # choose target DataFrame (do not use `df or self.rxns` because pandas.DataFrame is not boolean)
        df_in = self.rxns if df is None else df
        # operate on a copy to avoid surprising in-place side effects for caller
        df_work = df_in.copy()

        if df_work.empty:
            log.warn("No reactions found to plot.")
            return

        cf = self.checkframe
        outfolder = basename or self.basename+"_dir"
        if not os.path.exists(outfolder):
            os.makedirs(outfolder, exist_ok=True)

        for idx, rxn in df.iterrows():
            frame = rxn["frame"]
            timestep = rxn["timestep"]
            rxnID = rxn["rxnID"]
            rxnCount = rxn["rxnCount"]
            hash_before = rxn["rxn_hash_before"]
            hash_after = rxn["rxn_hash_after"]

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
            plt.suptitle(f"Timestep: {timestep} RxnType: {rxnID}\n {hash_before}:{hash_after}")
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
            pos_before = nx.spring_layout(Gbefore, iterations=75, pos=pos)
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
            pos_after = nx.spring_layout(Gafter, iterations=75, pos=pos)
            # plot data
            #nx.draw_networkx_edges(Gafter, edgelist=active_edges_after, alpha=0.6, width=5.0, edge_color="tab:red", pos=pos)
            nx.draw_networkx_edges(Gafter, edgelist=active_edges, alpha=0.6, width=5.0, edge_color="tab:red", pos=pos_after)
            nx.draw(Gafter, pos=pos_after, node_color=color, 
                    with_labels=True, labels=node_labels, font_size=6,
                    node_size=300, edge_color="black", width=bo)
            plt.tight_layout()
            fig = plt.gcf()
            fileout = f"{self.basename}_timestep{timestep}_rxnType{rxnID}_rxnCount{rxnCount}.png"
            f_out = os.path.join(outfolder, fileout)
            plt.savefig(f_out,dpi=200)
            plt.close(fig)

## analyze topology ##
# get degrees #
def get_degrees(df:pd.core.frame.DataFrame=None, target_atoms:tuple[int|str,...]=None ) -> list[list[int]]:
    df_work = df.copy()
    degrees = [None]*len(df_work)
    for idx, (df_idx, frame) in enumerate(df_work.iterrows()):
        nxg = frame["graph"]

        if target_atoms is None:
            nodes = list(nxg.nodes())
        else:
            nodes = list(node for node, node_data in nxg.nodes(data=True)
                            if node_data.get('type') in target_atoms 
                            or node_data.get('element') in target_atoms)

        degrees[idx] = list(d for n, d in nxg.degree(nodes))
    return degrees

    
def find_minimum_cycle_basis(df: pd.DataFrame = None, min_size: int = 7, max_block_size: int = None) -> list[list[int]]:
    """
    Computes the Minimum Cycle Basis (MCB) for each frame in the trajectory.
    This method identifies the Smallest Set of Smallest Rings (SSSR) by 
    decomposing the graph into biconnected components. 

    Args:
        df (pd.DataFrame, optional): Input DataFrame containing 'graph' column. 
            Defaults to self.backbone or self.frames.
        min_size (int): Minimum number of nodes for a cycle to be included.
        max_block_size (int, optional): Safety threshold. Blocks with more nodes 
            than this will be skipped to avoid O(n^3) complexity stalls.

    Returns:
        list: A nested list [frames][cycles][node_ids].

    Complexity Note:
    The MCB algorithm is O(m^3 * n). For dense networks that 'gel' into a 
    single large block, max_block_size is highly recommended to avoid 
    computational stalls.
    """
    # Select data source
    df_in = df
    all_frame_cycles = [None] * len(df_in)

    for idx, (df_idx, frame) in enumerate(df_in.iterrows()):
        nxg = frame["graph"]
        frame_basis = []

        # Use biconnected components to isolate cyclic parts of the network
        for block_nodes in nx.biconnected_components(nxg):
            block_len = len(block_nodes)
        
            if block_len < min_size:
                continue
        
            # Check safety limit for computational cost
            if max_block_size and block_len > max_block_size:
                log.warning(
                    f"Skipping large block ({block_len} nodes) in frame {idx}. "
                    f"Increase max_block_size if this analysis is required."
                )
                continue
            
            subgraph = nxg.subgraph(block_nodes)
        
            # MCB calculation (the heavy lifting)
            block_basis = nx.minimum_cycle_basis(subgraph)
        
            # Filter and store results
            frame_basis.extend([c for c in block_basis if len(c) >= min_size])
    
        all_frame_cycles[idx] = frame_basis
    
    return all_frame_cycles

