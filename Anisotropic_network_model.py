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

def ANM(filename):
    cuttoff = 12.0
    res_idx, ca_coordinate = read_ca_coordinates(filename)
    postall = ca_coordinate
    #print(postall)
    n2 = len(ca_coordinate)
    hessian = np.zeros((3*n2,3*n2),dtype=float)
    #print(netmat)
    for i in range(n2):
        for j in range(n2):
            if i==j:
                continue
            else:
                dis = ((ca_coordinate[i][0]-ca_coordinate[j][0]) ** 2 + (ca_coordinate[i][1]-ca_coordinate[j][1]) ** 2 + (ca_coordinate[i][2]-ca_coordinate[j][2]) ** 2) ** 0.5
                if dis <= cuttoff:
                    hessian[3 * i][3 * j] = -(postall[j][0] - postall[i][0]) * (postall[j][0] - postall[i][0]) / (dis * dis)
                    hessian[3 * i][3 * j + 1] = -(postall[j][0] - postall[i][0]) * (postall[j][1] - postall[i][1]) / (dis * dis)
                    hessian[3 * i][3 * j + 2] = -(postall[j][0] - postall[i][0]) * (postall[j][2] - postall[i][2]) / (dis * dis)
                    hessian[3 * i + 1][3 * j] = -(postall[j][1] - postall[i][1]) * (postall[j][0] - postall[i][0]) / (dis * dis)
                    hessian[3 * i + 1][3 * j + 1] = -(postall[j][1] - postall[i][1]) * (postall[j][1] - postall[i][1]) / (dis * dis)
                    hessian[3 * i + 1][3 * j + 2] = -(postall[j][1] - postall[i][1]) * (postall[j][2] - postall[i][2]) / (dis * dis)
                    hessian[3 * i + 2][3 * j] = -(postall[j][2] - postall[i][2]) * (postall[j][0] - postall[i][0]) / (dis * dis)
                    hessian[3 * i + 2][3 * j + 1] = -(postall[j][2] - postall[i][2]) * (postall[j][1] - postall[i][1]) / (dis * dis)
                    hessian[3 * i + 2][3 * j + 2] = -(postall[j][2] - postall[i][2]) * (postall[j][2] - postall[i][2]) / (dis * dis)
                    #print(dis)
                    #print( hessian[3 * i + 2][3 * j + 2])
                    hessian[3 * j][3 * i] = hessian[3 * i][3 * j]
                    hessian[3 * j + 1][3 * i] = hessian[3 * i][3 * j + 1]
                    hessian[3 * j + 2][3 * i] = hessian[3 * i][3 * j + 2]
                    hessian[3 * j][3 * i + 1] = hessian[3 * i + 1][3 * j]
                    hessian[3 * j + 1][3 * i + 1] = hessian[3 * i + 1][3 * j + 1]
                    hessian[3 * j + 2][3 * i + 1] = hessian[3 * i + 1][3 * j + 2]
                    hessian[3 * j][3 * i + 2] = hessian[3 * i + 2][3 * j]
                    hessian[3 * j + 1][3 * i + 2] = hessian[3 * i + 2][3 * j + 1]
                    hessian[3 * j + 2][3 * i + 2] = hessian[3 * i + 2][3 * j +2]

    for i in range(n2):
        for j in range(n2):
            if j!=i:
                dis = ((ca_coordinate[i][0] - ca_coordinate[j][0]) ** 2 + (
                            ca_coordinate[i][1] - ca_coordinate[j][1]) ** 2 + (
                                   ca_coordinate[i][2] - ca_coordinate[j][2]) ** 2) ** 0.5
                if dis<=cuttoff:
                    hessian[3 * i][3 * i] = hessian[3 * i][3 * i] + (postall[j][0] - postall[i][0]) * (postall[j][0] - postall[i][0]) / (dis * dis)
                    hessian[3 * i][3 * i + 1] = hessian[3 * i][3 * i + 1] + (postall[j][0] - postall[i][0]) * (postall[j][1] - postall[i][1]) / (dis * dis)
                    hessian[3 * i][3 * i + 2] = hessian[3 * i][3 * i + 2] + (postall[j][0] - postall[i][0]) * (postall[j][2] - postall[i][2]) / (dis * dis)
                    hessian[3 * i + 1][3 * i] = hessian[3 * i + 1][3 * i] + (postall[j][1] - postall[i][1]) * (postall[j][0] - postall[i][0]) / (dis * dis)
                    hessian[3 * i + 1][3 * i + 1] = hessian[3 * i + 1][3 * i + 1] + (postall[j][1] - postall[i][1]) * (postall[j][1] - postall[i][1]) / (dis * dis)
                    hessian[3 * i + 1][3 * i + 2] = hessian[3 * i + 1][3 * i + 2] + (postall[j][1] - postall[i][1]) * (postall[j][2] - postall[i][2]) / (dis * dis)
                    hessian[3 * i + 2][3 * i] = hessian[3 * i + 2][3 * i] + (postall[j][2] - postall[i][2]) * (postall[j][0] - postall[i][0]) / (dis * dis)
                    hessian[3 * i + 2][3 * i + 1] = hessian[3 * i + 2][3 * i + 1] + (postall[j][2] - postall[i][2]) * (postall[j][1] - postall[i][1]) / (dis * dis)
                    hessian[3 * i + 2][3 * i + 2] = hessian[3 * i + 2][3 * i + 2] + (postall[j][2] - postall[i][2]) * (postall[j][2] - postall[i][2]) / (dis * dis)



    #print(hessian)
    evalues, evectors = np.linalg.eig(hessian)
    #print(np.sort(evalues))
    ###sort the evalues and the corresponding evectors
    idx = evalues.argsort()#[::-1]
    evalues = evalues[idx]
    evectors = evectors[:,idx]
    print(evalues[1:30])
    return res_idx, n2, postall, evalues, evectors

