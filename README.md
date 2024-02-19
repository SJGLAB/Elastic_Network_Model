# Elastic_Network_Model
By Prof. Ji Guo Su
****
Repository for the application of elastic network mode, including Gaussian network model (GNM) and anisotropic network model (ANM), to extract the intrinsic dynamics, motion correlation, and entropy transfer encoded in the native structure of RSV prefusion protein trimer. 

****
- [Functions](#Functions)
- [Usage](#Usage)

## Functions
***Gaussian_nework_model.py***  constructs the Kirchhoff matrix and calculates the eigenvalues and eigenvectors.  
***Anisotropic_network_model.py***  constructs the Hessian matrix and calculates the eigenvalues and eigenvectors.  
***Flu_disp.py***  computes the residue fluctuations of the specified normal modes from GNM.  
***Cross_correlation.py***  computes the normalized cross-correlations of the fluctuations between different regions by the specified normal modes from GNM.  
***Mode_disp.py***  displays the residue motion of the specified normal modes from ANM.  
***Entropy_transfer.py***  is to analyze the entropy transfer in the structure of RSV pre-F trimer.  

## Usage
* To calculate the residue fluctuations using the first seven slowest normal modes from the GNM.
```
python3 Flu_disp.py 7
```
* To calculate the residue fluctuations using the first 30 slowest normal modes from the GNM.
```
python3 Flu_disp.py 30
```
* To calculate the normalized cross-correlation of the fluctuations between different regions using the first seven slowest normal modes from GNM.
```
python3 Cross_correlation.py 7
```
* To calculate the normalized cross-correlation of the fluctuations between different regions using the first 30 slowest normal modes from GNM.
```
python3 Cross_correlation.py 30
```
* To display the first slowest symmetric normal modes from the ANM.
```
python3 Mode_disp.py 8   #which generate a file named “modedisp_8.bild”, then a porcupine plot of the first slowest symmetric mode can be generated using UCSF Chimera software.
```
*  To display the second slowest symmetric normal modes from the ANM.
``` 
python3 Mode_disp.py 11   #which generate a file named “modedisp_11.bild”, then a porcupine plot of the second slowest symmetric mode can be generated using UCSF Chimera software.
```
* To calculate the entropy out from each residue with the parameter T/T_0 = 20.
```
python3 Entropy_transfer.py 20   #which generate two figures, and the second figure displays the results of the entropy out from each residue.
```
* To calculate the entropy out from each residue with the parameter T/T_0 = 10 (or 50).
```
python3 Entropy_transfer.py 10 (or 50)   #which generate two figures, and the second figure displays the results of the entropy out from each residue.
```
* To calculate the net entropy transfer between different regions with the parameter T/T_0 = 10 (or 50).
```
python3 Entropy_transfer.py 10 (or 50)   #which generate two figures, and the first figure displays the results of the net entropy transfer between different regions.

```
