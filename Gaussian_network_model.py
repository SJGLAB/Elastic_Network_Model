
import pdb
import numpy as np
import matplotlib.pyplot as plt

from Bio import PDB

#matplotlib.use('TkAgg')
def read_ca_coordinates(filename):
    parser = PDB.PDBParser()
    structure = parser.get_structure("protein", filename)
    ca_coordinate = []
    res_idx = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' not in residue: #or len(residue['CA']) == 0:
                    continue
                ca_coordinate.append(residue['CA'].coord)
                res_idx.append(residue.get_id()[1])

    print(len(ca_coordinate))
    return res_idx,ca_coordinate

def GNM(filename):
    cuttoff = 7.5
    res_idx, ca_coordinate = read_ca_coordinates(filename)
    n2 = len(ca_coordinate)
    netmat = np.zeros((n2,n2),dtype=int)
    #print(netmat)
    for i in range(n2):
        for j in range(n2):
            if i==j:
                continue
            else:
                dis = ((ca_coordinate[i][0]-ca_coordinate[j][0]) ** 2 + (ca_coordinate[i][1]-ca_coordinate[j][1]) ** 2 +  (ca_coordinate[i][2]-ca_coordinate[j][2]) ** 2) ** 0.5
                if dis <= cuttoff:
                    netmat[i][j] = -1
                    netmat[j][i] = -1

    for i in range(n2):
        netmat[i][i] = -1 * np.sum(netmat[i])
    #print(netmat)
    evalues, evectors = np.linalg.eig(netmat)
    #print(np.sort(evalues))
    ###sort the evalues and the corresponding evectors
    idx = evalues.argsort()#[::-1]
    evalues = evalues[idx]
    evectors = evectors[:,idx]
    print(evalues)
    return res_idx, n2, evalues, evectors

