from Gaussian_network_model import GNM
import numpy as np
import matplotlib.pyplot as plt
import sys

def normalized_cross_correlation(filename,mode_num):
    res_idx, n2, evalues, evectors = GNM(filename)
    cross_corr = np.zeros((n2,n2),dtype=float)
    norm_cross_corr = np.zeros((n2,n2),dtype=float)
    for i in range(n2):
        for j in range(n2):
            for k in range(mode_num):
                cross_corr[i][j] = cross_corr[i][j] + evectors[i][k+1] * evectors[j][k+1] / evalues[k+1]
    for i in range(n2):
        for j in range(n2):
            norm_cross_corr[i][j] = cross_corr[i][j] / np.sqrt(cross_corr[i][i] * cross_corr[j][j])
    res_num = 449
    norm_cross_corr_trueidx = np.zeros((res_idx[res_num-1]+1,res_idx[res_num-1]+1),dtype=float)
    for i in range(res_num):
        for j in range(res_num):
            norm_cross_corr_trueidx[res_idx[i]][res_idx[j]] = norm_cross_corr[i][j]
    plt.figure(num=1)
    plt.imshow(norm_cross_corr_trueidx.T, cmap='jet', origin='lower', interpolation=None)
    plt.colorbar()
    plt.show()


if __name__ == "__main__":
    filename = "4jhw.pdb"
    args = int(sys.argv[1])
    normalized_cross_correlation(filename,args)