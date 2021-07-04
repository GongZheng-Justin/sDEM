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

```
