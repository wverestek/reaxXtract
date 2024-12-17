#!/usr/bin/env python3
import sys, os.path, random
import networkx as nx
import matplotlib.pyplot as plt


startstep = 0
checkstep = 1
framestep = 1

rxn_bond_cutoff = 0
print_bonds = 5

type2color = {1:'grey',2:'antiquewhite',3:'antiquewhite',4:'red', 5:'red', 6:'red', 7:'red', 8:'red'}

# for reproducability of the plots
random.seed(a=23423)


# check number of arguments
if len(sys.argv)<2:
    print('ERROR: no file given')
    sys.exit(0)
    
# check filename
fname = sys.argv[1]
print('searching:',fname, end='')
if os.path.isfile(sys.argv[1]):
    print('\tfound')
else:
    print('\tnot found found')
    print('exiting')
    sys.exit(0)

# for output
outfolder = fname+'.dir'
if not os.path.exists(outfolder):
    os.makedirs(outfolder)

###################################
# get all neighbors of grade grad #
###################################
def getNeighsNgrade(G,nodes,grad):
    # make sure "nodes" is a list
    if isinstance(nodes, int): result = list([nodes])
    else: result = list(nodes)
    # find neighbors
    if grad>0:
        if isinstance(nodes, int):
            result = result + list(G[nodes])
        elif len(nodes)>=1:
            for n in nodes: result = result + list(G[n])
        else:
            print('error in neighbor search: empty nodes')
    # recursion
    if grad > 1:
        for n in result:
            result = result + getNeighsNgrade(G,n,grad-1)
    # flatten and remove duplicates
    if len(result)>1:
        result = list(set(result))
    return result


# f     filepointer
# fn    filename
# idx   running index
# ts    timestep
# pnum    number of particles
# mol   mol ID



### reading file
print('reading', fname)
f = open(fname,'r')
idx=-1
ts = []
pnum = []
nxg = []
line = f.readline()
while line:
    varline = line.strip().split()
    # Timestep
    if line.startswith('# Timestep'):
        idx = idx + 1
        nxg.append(nx.Graph())                      # array of networkx graphs
        ts.append(int(varline[-1]))                 # array of timesteps
        print('Reading timestep:',ts[idx],idx)
        line = f.readline()
        continue
    # Number of particles
    elif line.startswith('# Number of particles'):
        pnum.append(int(varline[-1]))
        line = f.readline()
        continue
    elif line.startswith('#') or  line.startswith('\n'):
        line = f.readline()
        continue
    elif len(varline)==0:
        line = f.readline()
        continue
    # lines with atom/bond info
    elif varline[0].isdigit():
        aidx = 0
        atomID = [0]*pnum[idx]
        atomType = [0]*pnum[idx]
        abo = [0.0]*pnum[idx]
        nlp = [0.0]*pnum[idx]
        q = [0.0]*pnum[idx]
        mol = [0]*pnum[idx]
        bonds = []
        
        while line:
            atomID[aidx] = int(varline[0])
            atomType[aidx] = int(varline[1])
            abo[aidx] = float(varline[-3])
            nlp[aidx] = float(varline[-2])
            q[aidx] = float(varline[-1])
            
            nb = int(varline[2])
            mol[aidx] = int(varline[3+nb])
            
            
            for tmp in range(nb):
                bonds.append( (atomID[aidx],int(varline[3+tmp]), float(varline[-4-tmp]) ) )

            line = f.readline()
                        # next Timestep
            if line.startswith('#') or line.startswith('\n'):
                break
            else:
                aidx = aidx + 1
                varline = line.strip().split()

        ###########################
        # fill Graph for timestep #
        ###########################
        # nodes
        tmp = [(a,{'AtomType':b,'abo':c,'nlp':d,'q':e}) for a,b,c,d,e in zip(atomID,atomType,abo,nlp,q) ]
        #print(tmp)
        nxg[idx].add_nodes_from(tmp)
        
        # edges 
        #tmp = [(a,b,{'bo':c,'weight':c}) for a,b,c in bonds]
        tmp = [(a,b,{'bo':c,'bos':c*3}) for a,b,c in bonds]
        #print(tmp)
        nxg[idx].add_edges_from(tmp)
        
        # conintue with next section '# Timestep' etc. or  continue below with processing

    else:
        print("You shouldn't be here...")
    
#########################
# finished reading file #
#########################
f.close()

##################
# find reactions #
##################
for idx in range(startstep, len(ts)-checkstep, framestep):
    print('Timestep',ts[idx],'-',ts[idx+checkstep],'reacting atoms:')
    #######################################################
    # comparing graphs, edges that are not in both graphs #
    #######################################################
    G1 = nx.difference(nxg[idx],nxg[idx+checkstep])
    G2 = nx.difference(nxg[idx+checkstep],nxg[idx])
    # atom IDs that are involved with changed bond connectivity
    reacting_atoms = list(set(list(sum(G1.edges, ())) + list(sum(G2.edges, ()))))
    del G1,G2
    
    # no reaction found
    if len(reacting_atoms) < 1: continue

    #####################################
    # include reactions rxn_bond_cutoff away? #
    #####################################
    if rxn_bond_cutoff > 0:
        env_atoms1 = set(getNeighsNgrade(nxg[idx],reacting_atoms,rxn_bond_cutoff))
        env_atoms2 = set(getNeighsNgrade(nxg[idx+checkstep],reacting_atoms,rxn_bond_cutoff))
        reacting_atoms = list(env_atoms1.union(env_atoms2))
        del env_atoms1, env_atoms2

    #########################
    # group by connectivity #
    #########################
    # subgraph before reaction
    Gchange = nxg[idx].subgraph(reacting_atoms).copy()
    # add edges after reaction
    Gchange.add_edges_from(nxg[idx+checkstep].subgraph(reacting_atoms).edges)
    # extract groups of connected atoms
    reaction_sets = [tmp for tmp in list(nx.connected_components(Gchange))]
    del Gchange
        
    # printing
    for tmp in reaction_sets:
        print(tmp)
    
    for i,rset in enumerate(reaction_sets):
        # get reacting edges
        G1 = nxg[idx].subgraph(rset)
        G2 = nxg[idx+checkstep].subgraph(rset)
        active_edges1 = nx.difference(G2,G1).edges()
        active_edges2 = nx.difference(G1,G2).edges()
        
        # reacting atoms + print environament atoms
        print_atoms = set(getNeighsNgrade(nxg[idx+checkstep],list(rset),print_bonds))
        Gprint =  nxg[idx].subgraph(print_atoms).copy()
        # positions with spring model
        Gprint.add_edges_from(active_edges1)
        pos = nx.spring_layout(Gprint,iterations=5000, dim=3, seed=random.randrange(0,5000))
        pos = {k:v[0:2] for k,v in pos.items()}
        #pos = nx.spring_layout(Gprint, iterations=5000, scale=1.0, seed=random.randrange(0,5000))
        #Gprint.add_edges_from(active_edges1)
        pos = nx.spring_layout(Gprint, iterations=5000, scale=1.0, pos=pos)
        #pos = nx.spring_layout(Gprint, iterations=1000, scale=1.0, pos=pos)
        # node color
        node2type = nx.get_node_attributes(Gprint,'AtomType')
        color = [type2color[node2type[i]] for i in Gprint.nodes]
        del Gprint

        # before reaction
        plt.clf()
        Gprint =  nxg[idx].subgraph(print_atoms)
        pos = nx.spring_layout(Gprint, iterations=5, pos=pos)
        bo = [tmp**2 for tmp in list(nx.get_edge_attributes(Gprint,'bo').values())]
        nx.draw_networkx_edges(Gprint, edgelist=active_edges2, width=5.0, edge_color='tab:red',pos=pos)
        nx.draw(Gprint, pos=pos, with_labels=True, node_color=color, edge_color='black', width=bo)
        fout = os.path.join(outfolder,fname+'_'+str(ts[idx+checkstep])+'_R'+str(i)+'_0.png')
        plt.savefig(fout)

        # after reaction
        plt.clf()
        Gprint =  nxg[idx+checkstep].subgraph(print_atoms)
        pos = nx.spring_layout(Gprint, iterations=5, pos=pos)
        bo = [tmp**2 for tmp in list(nx.get_edge_attributes(Gprint,'bo').values())]
        nx.draw_networkx_edges(Gprint, edgelist=active_edges1, width=5.0, edge_color='tab:red',pos=pos)
        nx.draw(Gprint, pos=pos, with_labels=True, node_color=color, edge_color='black', width=bo)
        fout = os.path.join(outfolder,fname+'_'+str(ts[idx+checkstep])+'_R'+str(i)+'_1.png')
        plt.savefig(fout)

    #del reaction_sets
    
