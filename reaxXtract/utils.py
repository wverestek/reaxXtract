import os
import warnings
import networkx as nx

###############################
# --- Static dictionaries --- #
###############################
# atomic ordinal number => element symbol
ON2ELEM: dict[str, str] = { 
    0: "X",
    1: "H" ,  2: "He",
    3: "Li",  4: "Be",  5: "B" ,  6: "C" ,  7: "N" ,  8: "O" ,  9: "F" , 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P" , 16: "S" , 17: "Cl", 18: "Ar",
    19: "K" , 20: "Ca", 21: "Sc", 22: "Ti", 23: "V" , 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 
    28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr"
    }

ON2MASS: dict[int, float] = {
     0: 12.011,
     1:  1.008,  2:  4.0026,   
     3:  6.940,  4:  9.0122,  5: 10.810,  6: 12.011,  7: 14.0070,  8: 15.9990,  9: 18.998, 10: 20.180,
    11: 22.990, 12: 24.3050, 13: 26.982, 14: 28.085, 15: 30.9740, 16: 32.0600, 17: 35.450, 18: 39.948, 
    19: 39.098, 20: 40.0780, 21: 44.956, 22: 47.867, 23: 50.9415, 24: 51.9961, 25: 54.938, 26: 55.845, 27: 58.933, 
    28: 58.6934, 29: 63.546, 30: 65.380, 31: 69.723, 32: 72.6300, 33: 74.9216, 34: 78.971, 35: 79.904, 36: 83.798
    }

ON2RAD: dict[int, float] = {
     0: 170.0,
     1: 110,  2: 140,   
     3: 182,  4: 153,  5: 192,  6: 170,  7: 155,  8: 152,  9: 147, 10: 154,
    11: 227, 12: 173, 13: 184, 14: 210, 15: 180, 16: 180, 17: 175, 18: 188, 
    19: 275, 20: 231, 21: 200, 22: 200, 23: 200, 24: 200, 25: 200, 26: 200, 27: 200, 
    28: 163, 29: 140, 30: 139, 31: 187, 32: 211, 33: 185, 34: 190, 35: 185, 36: 202
    } # in pm, 21:Sc to 27:Co are dummy values


# default color fallback
DEFAULT_COLOR = "#E0E0E0"

# atomic ordinal number => RGB color (hex)
# Jmol color scheme: https://jmol.sourceforge.net/jscolors/
# antiquewhite FAEBD7 instead of white FFFFFF for Hydrogen
ON2HEX: dict[int, str] = { 
    0: "E0E0E0",
    1: "#FAEBD7",  2: "#D9FFFF",  3: "#CC80FF",  4: "#C2FF00",  5: "#FFB5B5",  6: "#909090",
    7: "#3050F8",  8: "#FF0D0D",  9: "#90E050", 10: "#B3E3F5", 11: "#AB5CF2", 12: "#8AFF00",
    13: "#8FAF00", 14: "#F0C8A0", 15: "#FF8000", 16: "#FFFF30", 17: "#1FF01F", 18: "#80D1E3",
    19: "#8F40D4", 20: "#3DFF00", 21: "#E6E6E6", 22: "#BFC2C7", 23: "#A6A6AB", 24: "#8A99C7", 25: "#9C7AC7", 26: "#E06633", 27: "#F090A0",
    28: "#50D050", 29: "#C88033", 30: "#7D80B0", 31: "#C28F8F", 32: "#668F8F", 33: "#BD80E3", 34: "#FFA100", 35: "#A62929", 36: "#5CB8D1"
    }


# element symbol => RGB color (hex)
ELEM2HEX:dict[str, str] = {ON2ELEM[i]:ON2HEX[i] for i in ON2ELEM}




#####################
# --- Utilities --- #
#####################

def k_nearest_neighs(G: nx.Graph, start: set[int], k: int) -> set[int]:
    """
    Return all nodes within k bonds of a given start set in a networkx Graph G.
    """
    neighs = set(start)
    for l in range(k):
        #neighs.update( set((nbr for n in neighs for nbr in G[n])) )
        neighs |=  {nbr for n in neighs for nbr in G[n]}
    return neighs


#def convert_str2dict(mystr: str) -> dict:
#    if len(mystr) > 0 and ":" in mystr:
#        return {k: v for k, v in [[int(i) for i in tmp.split(":")] for tmp in mystr.strip().split(",")]}
#    else:
#        return dict()
def convert_str2dict(mystr: str) -> dict[int, int]:
    """
    Convert a mapping string like "1:6,2:8,3:1" into a dict {1:6, 2:8, 3:1}.
    - Whitespace-tolerant
    - Issues warnings for malformed entries
    - Returns {} for empty or invalid input
    """
    result = {}
    if not mystr:
        return result
    
    for pair in mystr.split(","):
        pair = pair.strip()
        #print (f"Processing pair: '{pair}'")
        if not pair or ":" not in pair:
            if pair:
                warnings.warn(
                    f"Ignoring malformed mapping entry (missing ':'): '{pair}'",
                    UserWarning
                )
            continue
        try:
            #k, v = (int(x.strip(),int(y.strip())) for x,y in pair.split(":", 1))
            #print(pair.split(":", 1))
            k, v = (x for x in pair.split(":", 1))
            result[int(k)] = v
        except ValueError:
            warnings.warn(
                f"Ignoring invalid mapping entry (non-integer): '{pair}'",
                UserWarning
            )
    return result



