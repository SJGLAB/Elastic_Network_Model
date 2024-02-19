from Gaussian_network_model import GNM
import numpy as np
import matplotlib.pyplot as plt
import sys

def fluctuation_dis(filename,mode_num):
    res_idx, n2, evalues, evectors = GNM(filename)
    flu = np.zeros((n2,1),dtype=float)
    for i in range(n2):
        for j in range(mode_num):
            flu[i] = flu[i] + evectors[i][j+1] * evectors[i][j+1] / evalues[j+1]
    plt.plot(res_idx[0:449],flu[0:449],color='k')
    plt.plot(res_idx[72:151], flu[72:151], color='b')
    plt.plot(res_idx[47:71], flu[47:71], color='r')
    plt.plot(res_idx[151:176], flu[151:176], color='r')
    #plt.legend()
    plt.show()

if __name__ == "__main__":
    filename = "4jhw.pdb"
    args = int(sys.argv[1])
    fluctuation_dis(filename,args)