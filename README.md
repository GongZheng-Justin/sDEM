# sDEM
&emsp;**sDEM** (simple Discrete Element Method) is a simple parallel DEM code with following features:

* MPI parallelization by means of 2D Domain Decomposition
* Have the ability to handle particles with different diameters
* Both linear and non-linear contact force models
* NBS-Munjiza (Non Binary Search - Munjiza) contact search algorithm

&emsp;Following picture shows the process of particles settling in a sand box.
![](doc/sandbox.gif)


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
&emsp;If the compiling processes successfully, the executable file(s) `dem`/`channel4th` will be appeared in the folder `sDEM-master/`.

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
## To do list :muscle:
* Adding a **non-spherical Particles** module, using **Super-ellipsoids** and/or **Multi-Sphere** Method
* Hybrid MPI/OpenMP parallelization and GPU acceleration  

## Acknowledgements :clap:
  I would particularly thank [Dr. Norouzi](https://www.researchgate.net/profile/Hamid-Norourzi) from University of Tehran, and his book **_Coupled CFD‐DEM Modeling: Formulation, Implementation and Applimation to Multiphase Flows_**, besides the [attached DEM code](https://www.wiley.com//legacy/wileychi/norouzi/form.html?type=SupplementaryMaterial). 

## Contact and Feedback :email:
&emsp;If you have any question, or want to contribute to the code, please don't hesitate to contact me: Zheng Gong (gongzheng_justin@outlook.com)
