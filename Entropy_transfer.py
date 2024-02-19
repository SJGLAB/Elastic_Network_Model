from Gaussian_network_model import GNM
import numpy as np
import math
import matplotlib.pyplot as plt
import sys

def entropytrans (filename,t_t0):
    res_idx, n2, evalues, evectors = GNM(filename)
    T = np.zeros((n2,n2),dtype=float)

    A_sumk = np.zeros((n2,n2),dtype=float)
    A_sumkexp = np.zeros((n2,n2),dtype=float)
    k_c = n2
    for i in range(n2):
        print(i)
        for j in range(n2):
            for k in range(1,k_c):
                A_sumk[i][j] = A_sumk[i][j] + evectors[i][k] * evectors[j][k] / evalues[k]
                A_sumkexp[i][j] = A_sumkexp[i][j] + (evectors[i][k] * evectors[j][k] / evalues[k]) * math.exp(-1 * evalues[k] * t_t0)

    for i in range(n2):
        for j in range(n2):
            if i != j:
                Tvalue = 0.5 * math.log(A_sumk[j][j] * A_sumk[j][j] - A_sumkexp[j][j] * A_sumkexp[j][j]) - 0.5 * math.log(A_sumk[i][i] * A_sumk[j][j] * A_sumk[j][j] + 2.0 * A_sumk[i][j] * A_sumkexp[j][j] * A_sumkexp[i][j] - (A_sumkexp[i][j] * A_sumkexp[i][j] + A_sumk[i][j] * A_sumk[i][j]) * A_sumk[j][j] - A_sumkexp[j][j] * A_sumkexp[j][j] * A_sumk[i][i]) - 0.5 * math.log(A_sumk[j][j]) + 0.5 * math.log(A_sumk[i][i] * A_sumk[j][j] - A_sumk[i][j] * A_sumk[i][j])
                print(Tvalue)
                T[i][j] = Tvalue

    delt_T = np.zeros((n2,n2),dtype=float)
    for i in range(n2):
        for j in range(n2):
            delt_T[i][j] = T[i][j]-T[j][i]

    nnc = 449
    delt_T_trueindex = np.zeros((res_idx[nnc-1]+1,res_idx[nnc-1]+1),dtype=float)
    for i in range(nnc):
        for j in range(nnc):
            delt_T_trueindex[res_idx[i]][res_idx[j]] = delt_T[i][j]

    plt.figure(num=1)
    plt.imshow(delt_T_trueindex.T, cmap='jet', origin='lower', interpolation=None)
    plt.colorbar()
    plt.show()

    TE = np.zeros((n2,1),dtype=float)
    for i in range(n2):
        TE[i] = np.sum(T[i])

    plt.figure(num=2)
    plt.plot(res_idx[0:449], TE[0:449], color='k')
    plt.plot(res_idx[72:151], TE[72:151], color='b')
    plt.plot(res_idx[47:71], TE[47:71], color='r')
    plt.plot(res_idx[151:176], TE[151:176], color='r')
    plt.show()


if __name__ == "__main__":
    filename = "4jhw.pdb"
    args = float(sys.argv[1])
    entropytrans(filename, args)



