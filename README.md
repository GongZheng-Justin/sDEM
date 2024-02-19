# sDEM
&emsp;**sDEM** (simple Discrete Element Method) is a simple parallel DEM code with following features:

* MPI parallelization by means of 2D Domain Decomposition
* Have the ability to handle particles with different diameters
* Both linear and non-linear contact force models
* NBS-Munjiza (Non Binary Search - Munjiza) contact search algorithm

## Notice :newspaper: 
&emsp;**sDEM** has been integrated into the CFD-DEM sovler [**CP3d**](https://github.com/GongZheng-Justin/CP3d) as a sub-solver.

## Installation :briefcase:
&emsp;Present solver has the following two prerequisities:

* MPI
* GFortran/Intel Fortran (Supporting Fortran 2003 or higher version)

&emsp;You can compile the code as follows:
```
1. chmod a+x ./mymake.sh
2. ./mymake.sh
3. choose the correct compiler you use, just following guidances printed in the terminal
```
&emsp;If the compiling processes successfully, the executable file `dem` will be appeared in the current folder (e.g. `sDEM-master`).

## Usage :book:
&emsp;After compiling the code successfully, you can run the executable file like that:
```
mpirun -n [np] ./dem [inputFile]
```
&emsp;Here:
* `np` denotes the number of processors you use
* `inputFile` is the name string for the input parameter file  

&emsp;For instance, if you want to run the particles settling in a sand box, you can type the following words in your terminal:
```
mpirun -n 4 ./dem ./Input/SandBox_DEM.prm
```

## Postprocessing :space_invader:
&emsp;`Paraview` is used as the postprocessing tool as follows:
```
1. cd DEM/Results
2. paraview &
3. use paraview to open the interface XDMF file (e.g. PartVisuForSandBox.xmf)
```

## To do list :muscle:
* Adding a **non-spherical Particles** module, using **Super-ellipsoids** and/or **Multi-Sphere** Method
* Hybrid MPI/OpenMP parallelization and GPU acceleration  

## Acknowledgements :clap:
&emsp;I would particularly thank [Dr. Norouzi](https://www.researchgate.net/profile/Hamid-Norourzi) from University of Tehran, for his continuous help (from superficial Fortran language explanation to the underlying DEM details), and his book **_Coupled CFD‐DEM Modeling: Formulation, Implementation and Applimation to Multiphase Flows_**, besides the attached code [cemfDEM](https://github.com/hamidrezanorouzi/cemfDEM). 

&emsp;To some extent, the present DEM solver can be regarded as the MPI version of the original serial [cemfDEM](https://github.com/hamidrezanorouzi/cemfDEM).

## Contact and Feedback :email:
&emsp;If you have any question, or want to contribute to the code, please don't hesitate to contact me: Zheng Gong (gongzheng_justin@outlook.com)

##
Following picture shows the process of particles settling in a sand box.  
<img src="./doc/sandbox.gif" width="30%" height="30%" div align=center />

## Star History
[![Star History Chart](https://api.star-history.com/svg?repos=GongZheng-Justin/sDEM&type=Date)](https://star-history.com/#GongZheng-Justin/sDEM&Date)

<img src="[![Star History Chart](https://api.star-history.com/svg?repos=GongZheng-Justin/sDEM&type=Date)](https://star-history.com/#GongZheng-Justin/sDEM&Date)" width="30%" height="30%" div align=center />
