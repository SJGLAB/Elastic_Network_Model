# Elastic_Network_Model
By Prof. Ji Guo Su
****
## Overview
Repository for the application of elastic network mode, including Gaussian network model (GNM) and anisotropic network model (ANM), to extract the intrinsic dynamics, motion correlation, and entropy transfer encoded in the native structure of RSV prefusion protein trimer. 
****
- [Functions](#Functions)
- [System_Requirements](#System_Requirements)
- [Installation_Guide](#Installation_Guide)
- [Usage_and_Demo](#Usage_and_Demo)
- [Instructions_for_use](#Instructions_for_use)

## Functions
***Gaussian_nework_model.py***  constructs the Kirchhoff matrix and calculates the eigenvalues and eigenvectors.  
***Anisotropic_network_model.py***  constructs the Hessian matrix and calculates the eigenvalues and eigenvectors.  
***Flu_disp.py***  computes the residue fluctuations of the specified normal modes from GNM.  
***Cross_correlation.py***  computes the normalized cross-correlations of the fluctuations between different regions by the specified normal modes from GNM.  
***Mode_disp.py***  displays the residue motion of the specified normal modes from ANM.  
***Entropy_transfer.py***  is to analyze the entropy transfer in the structure of RSV pre-F trimer.  

## System_Requirements
### **1. Hardware requirements**  

Only a standard computer with enough RAM to support the in-memory operations is required.  

### **2. Software requirements**  

#### OS requirements  
The codes are tested on the following OSes:   
- Linux x64  
- Windows 10 x64

And the following x86_64 version of Python:  
- Python 3.10
  
#### Python dependencies  
- pdb  
- numpy   
- matplotlib
- Bio 
- sys  
- math

## Installation_Guide
### Download the codes
```
git clone https://github.com/SJGLAB/Elastic_Network_Model.git
```
### Prepare the environment
We recommand you to use [Anaconda](https://www.anaconda.com/) to prepare the environments.  
You can create the environment via  
```
conda create -n Elastic_net python=3.10
conda activate Elastic_net
cd Elastic_Network_Model
pip install -r requirements.txt
```
The installation takes about 10 mins. 

## Usage_and_Demo
* To calculate the residue fluctuations using the first seven slowest normal modes from the GNM.
```
python3 Flu_disp.py 7   #The output is the residue fluctuation plot. This takes about several mins on a normal desktop computer.
```
* To calculate the residue fluctuations using the first 30 slowest normal modes from the GNM.
```
python3 Flu_disp.py 30   #The output is the residue fluctuation plot. This takes about several mins on a normal desktop computer.
```
* To calculate the normalized cross-correlation of the fluctuations between different regions using the first seven slowest normal modes from GNM.
```
python3 Cross_correlation.py 7   # The output is the normalized cross-correlation figure. This takes about several mins on a normal desktop computer.
```
* To calculate the normalized cross-correlation of the fluctuations between different regions using the first 30 slowest normal modes from GNM.
```
python3 Cross_correlation.py 30   # The output is the normalized cross-correlation figure. This takes about several mins on a normal desktop computer.
```
* To display the first slowest symmetric normal modes from the ANM.
```
python3 Mode_disp.py 8   #which generate a file named “modedisp_8.bild”, then a porcupine plot of the first slowest symmetric mode can be generated using UCSF Chimera software. This takes about several mins on a normal desktop computer.
```
* To display the second slowest symmetric normal modes from the ANM.
``` 
python3 Mode_disp.py 11   #which generate a file named “modedisp_11.bild”, then a porcupine plot of the second slowest symmetric mode can be generated using UCSF Chimera software. This takes about several mins on a normal desktop computer.
```
* To calculate the entropy out from each residue with the parameter T/T_0 = 20.
```
python3 Entropy_transfer.py 20   #which generate two figures, and the second figure displays the results of the entropy out from each residue. This takes about an hour on a normal desktop computer.
```
* To calculate the net entropy transfer between different regions with the parameter T/T_0 = 20.
```
python3 Entropy_transfer.py 20   #which generate two figures, and the first figure displays the results of the net entropy transfer between different regions.
```
* To calculate the entropy out from each residue with the parameter T/T_0 = 10 (or 50).
```
python3 Entropy_transfer.py 10 (or 50)   #which generate two figures, and the second figure displays the results of the entropy out from each residue. This takes about an hour on a normal desktop computer.
```
* To calculate the net entropy transfer between different regions with the parameter T/T_0 = 10 (or 50).
```
python3 Entropy_transfer.py 10 (or 50)   #which generate two figures, and the first figure displays the results of the net entropy transfer between different regions. This takes about an hour on a normal desktop computer.
```

## Instructions_for_use
The above **Usage and Demo** reproduces all the quantitative results in the manuscript.

## Licence
The code is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-nc-sa/4.0/).

## If you use this code, please cite the paper
Yu Liang, Shuai Shao, Xin Yu Li, Zi Xin Zhao, Ning Liu, Zhao Ming Liu, Fu Jie Shen, Hao Zhang, Jun Wei Hou, Xue Feng Zhang, Yu Qin Jin, Li Fang Du, Xin Li, Jing Zhang, Ji Guo Su & Qi Ming Li. 2024.

