# Energy-Density-Class
Energy Density Class
This class describes the subsequent evolution of the energy density (ED) and flow velocity (FV) of a number of partons after a heavy particle collision. One object of the class describes one system of partons and the information on the system's evolution is output into .txt files. The ED and FV are calculated on the (x, y) grid transverse to the particle collision. Each position on the grid corresponds to an ED value and a FV value at a given tau.
## Parameters ##
A given object (system of partons) requires the following parameters before the energy density and flow velocity values can be calculated:

**Number of partons and their positions**

```int partNum``` is the integer number of partons in the system and ```vector<pair<double, double>> partPos``` lists the positions of each parton. To add a parton, the user adds a new pair of x and y values to ```partPos```. ```partNum``` is updated by calculating the new size of the vector.

**Sigma**

The standard deviation in the Gaussian distribution part of the energy-momentum tensor

**τ<sub>0</sub>, τ<sub>f</sub>, τ<sub>step</sub>**

The variable τ, the proper time, is the analog of time for the evolution of this system. ```double tau0``` is the initial τ from  which the system's evolution begins. ```double tauFinal``` is the final time at which the system's evolution ends, and ```double tauStep``` is the time at which the evolution progresses. The system will be examined at (```tauFinal``` - ```tau0``` - 1) / ```tauStep``` points in time.

**x range, x scale, y range, y scale**

The size and complexity of the transverse grid on which the system is analyzed needs to be specified. The range of values in the x and y directions respectively are ```double xRange``` multiplied by ```double xScale``` and ```double yRange``` multiplied by ```double yScale```. The The transverse grid always starts at -1 * ```xRange * xScale``` and ends at ```xRange * xScale``` with steps in position differing by ```xScale``` for the x coordinates. The same is true for the y coordinates.
For instance, suppose ```xRange``` and ```yRange``` are both set to 12 and ```xScale``` and ```yScale``` are both set to 1. The points at the far right start at 12 and the points at the far left end at -12: (12, y), (11, y), (10, y)... (-11, y), (-12, y). This particular example makes for ((xRange * 2) + 1) * ((yRange * 2) + 1)  = ((12 * 2) + 1) * ((12 * 2) + 1) = 625 different points on the grid to be analyzed.
## Calculations ##
To calculate the numerous ED and FV values at every τ value and position, the eigenvectors and eigenvalues of the energy-momentum tensor (EMT) need to computed.

##Energy-momentum tensor calculator function: ```Eigen::MatrixXd getEMTensor(double x, double y, double tau)```##

This function computes the EMT at an x, y, and tau value passed to it. The returned matrix is the sum of the EMTs from each parton present in the system since the EMT function's parameters include the x and y values of a single particle, not multiple particles.
##Energy density calculator function: ```double energyDensity::getEnergyDensity(double x, double y, double tau)```##

After the EMT was calculated, the ED value can be found by retrieving the first eigenvalue from the EMT. This function returns the ED value at a given x, y, and τ.

##The flow velocity functions##

The FV in the x-direction is equal to the second component of the eigenvector and the FV in the y-direction is equal to the third component of the eigenvector. ```pair <double, double> energyDensity::getFlowVelocityXY(double x, double y, double tau)``` retrieves these values and returns both values at the same time as a pair. ```double energyDensity::getFlowVelocityX(double x, double y, double tau)``` calls ```getFlowVelocityXY(double x, double y, double tau)``` and returns the first value in the pair returned from the XY function. The same function, ```getFlowVelocityY(double x, double y, double tau)``` and process exists for the FV's y component.

##Energy density output function at single time: ```void energyDensity::EDGrid(double tau)```##

This function calls the ```getEnergyDensity(double x, double y, double tau)``` function for each point on the grid and outputs a .txt file with the ED value for each point, where the x and y values of each point are printed followed by their ED value. Remember that this function prints the ED for every point on the grid at a single point in time.

##Energy density output function at single position at all times: ```void energyDensity::EDEvolution(double x, double y)```##

This function calls the ```getEnergyDensity(double x, double y, double tau)``` function for a ***single point on the grid***, once for each point in time between ```tau0``` and ```tauFinal``` with a step of ```tauStep```.

##Energy density output function at all positions at all times: ```void energyDensity::EDEvolution()```##

This function calls the ```getEnergyDensity(double x, double y, double tau)``` function for **all points on the grid**, once for each point in time between ```tau0``` and ```tauFinal``` with a step of ```tauStep```

##Flow velocity evolution##
