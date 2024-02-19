from Anisotropic_network_model import ANM
import numpy as np
import matplotlib.pyplot as plt
import sys

def mode_display(filename,mode_idx):
    res_idx, n2, postall, evalues, evectors = ANM(filename)
    postall_add = np.zeros((n2,3),dtype=float)
    amplitude = 25.0
    fn = 'modedisp_' + str(mode_idx) + '.bild'
    file = open(fn, "w")
    for i in range(n2):
        postall_add[i][0] = postall[i][0] - (amplitude / evalues[mode_idx]) * evectors[3 * i][mode_idx]
        postall_add[i][1] = postall[i][1] - (amplitude / evalues[mode_idx])* evectors[3 * i + 1][mode_idx]
        postall_add[i][2] = postall[i][2] - (amplitude / evalues[mode_idx]) * evectors[3 * i + 2][mode_idx]

        content = '.arrow ' + str(postall[i][0]) + ' ' + str(postall[i][1]) + ' ' + str(postall[i][2]) + ' ' + str(postall_add[i][0]) + ' ' + str(postall_add[i][1]) + ' ' + str(postall_add[i][2]) + ' ' + ' 0.1 0.4 0.75' + '\n'
        file.write(content)
    file.close()

if __name__ == "__main__":
    filename = "4jhw.pdb"
    args = int(sys.argv[1])
    mode_display(filename,args)