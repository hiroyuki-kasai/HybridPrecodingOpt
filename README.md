# HybridPrecodingOpt : Optimization algorithms for hybrid precoding in millimeter wave (mmWave) MIMO systems
----------

Authors: [Hiroyuki Kasai](http://kasai.kasailab.com/)

Last page update: July 02, 2018

Latest version: 1.0.0 (see Release notes for more info)

<br />

Introduction
----------
Hybrid [precoding](https://en.wikipedia.org/wiki/Precoding)(-[beamforming](https://en.wikipedia.org/wiki/Beamforming)) is the most promising approach to reduce high hardware costs and high power consumptions in large-scale millimeter wave (mmWave) [MIMO](https://en.wikipedia.org/wiki/MIMO) systems. Hybrid precoding combines large-dimensional analog precoding (or beamforming) via phase shifters with lower-dimensional digital baseband precoding.

The maximization of spectral efficiency approximately boils down to a minimization problem of the Euclidean distance between the fully digital precoder and the hybrid precoder, 
This problem is further formulated as a [matrix factorization](https://it.wikipedia.org/wiki/Matrix_factorization) problem of he fully digital precoder with a product of the digital baseband precoder matrix and the analog RF precoder (or beamforming) matrix. 
The noteworthy point is that the phase shifters impose an additional element-wise unit modulus constraints on the analog RF precoder matrix. 

This package provides the codes of the proposed optimization algorithms for hybrid precoding. This code includes existing state-of-the arts algorithms, too. Most of the codes of this package come from the [brilliant project](https://github.com/yuxianghao/Alternating-minimization-algorithms-for-hybrid-precoding-in-millimeter-wave-MIMO-systems).  


<br />

Document
----------
The document can be found below;

- H. Kasai, "Fast optimization algorithm on complex oblique manifold for hybrid precoding in Millimeter Wave MIMO systems," arXiv paper (coming soon!).

<br />

Algorithms
----------

- **Proposed**
    - H. Kasai, "Fast optimization algorithm on complex oblique manifold for hybrid precoding in Millimeter Wave MIMO systems," arXiv paper (coming soon!).


- **MO-AltMin**
    - X. Yu, J.-C. Shen, J. Zhang, and K. B. Letaief, "[Alternating minimization algorithms for hybrid precoding in millimeter wave MIMO systems](https://ieeexplore.ieee.org/document/7397861/)," IEEE Journal on Selected Areas in Communications, vol. 10, no. 3, pp. 485-500, 2016.

- **PE-AltMin**
    - X. Yu, J.-C. Shen, J. Zhang, and K. B. Letaief, "[Alternating minimization algorithms for hybrid precoding in millimeter wave MIMO systems](https://ieeexplore.ieee.org/document/7397861/)," IEEE Journal on Selected Areas in Communications, vol. 10, no. 3, pp. 485-500, 2016.

- **OMP-based**
    - O. E. Ayach, S. Rajagopal, S. Abu-Surra, Z. Pi, and R. W. Heath, "[Spatially sparse precoding in millimeter wave MIMO systems](https://ieeexplore.ieee.org/document/6717211/)," IEEE Transations on Wireless Communications, vol. 13, no. 3, pp. 1499-1513, 2014.

<br />


Folders and files
---------
<pre>
./                      - Top directory.
./README.md             - This readme file.
./run_me_first.m        - The scipt that you need to run first.
./demo.m                - Demonstration script. 
./comp_OFDM.m           - Simulation script for OFDM. 
|proposed/              - Contains the files for the proposed algorithms.
|AltMinAlg/             - Contains the files for the existing algorithms.
|cvx/                   - Folder for CVX project (Please downlod by yoursel!).
</pre>

<br />  

First to do
----------------------------
Run `run_me_first` for path configurations. 
```Matlab
%% First run the setup script
run_me_first; 
```

<br />

Demo example for OFDM system
----------------------------

Just execute `demo` for the first demonstration of this package. 

<div align="center"><img src="http://www.kasailab.com/public/github/HybridPrecodingOpt/images/demo.png" width="500"></div>


<br />

Simulation for narrow band system
----------------------------

Execute `channel_realization` in `AltMinAlg/datasets/` folder to generate datasets such as `Ns=3.mat'. The script outputs the datasete file in there. You can load this generated dataset file in the script. 

<br />


More results for OFDM system
----------------------------

- Comparison with the MO-AltMin algorithm. 

<div align="center"><img src="http://www.kasailab.com/public/github/HybridPrecodingOpt/images/Comp_MOAltMin.png" width="900"></div>

<br />
<br />

- Spectral efficiency for SNRs. 

<div align="center"><img src="http://www.kasailab.com/public/github/HybridPrecodingOpt/images/OFDB_SNR.png" width="800"></div>

<br />
<br />

- Spectral efficiency and processing time under different number of data streams, and antennass in transmitter/receiver.

<div align="center"><img src="http://www.kasailab.com/public/github/HybridPrecodingOpt/images/OFDM_NRF.png" width="800"></div>

<br />
<br />

Notes
-------
- Most of the codes of the package come from the codes in [AltMinAlg](https://github.com/yuxianghao/Alternating-minimization-algorithms-for-hybrid-precoding-in-millimeter-wave-MIMO-systems).  
- The project uses the MATLAB toolbox [Manopt](https://www.manopt.org/).
- [CVX package](http://cvxr.com/cvx/) is required for the narrrowband simulations. Please download the packege from [here](http://cvxr.com/cvx/download/) into `cvx` folder. 

<br />


Problems or questions
---------------------
If you have any problems or questions, please contact the author: [Hiroyuki Kasai](http://kasai.kasailab.com/) (email: kasai **at** is **dot** uec **dot** ac **dot** jp)

<br />

Release Notes
--------------
* Version 1.0.0 (July 01, 2018)
