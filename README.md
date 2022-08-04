# Energy-Density-Class
Energy Density Class
This class describes the subsequent evolution of the energy density (ED) and flow velocity (FV) of a number of partons after a heavy particle collision. One object of the class describes one system of partons and the information on the system's evolution is output into .txt files. The ED and FV are calculated on the (x, y) grid transverse to the particle collision. Each position on the grid corresponds to an ED value and a FV value at a given tau.

# Parameters 
A given object (system of partons) requires the following parameters before the energy density and flow velocity values can be calculated:
- Number of partons and their positions
````int partNum```` is the integer number of partons in the system and ````vector<pair<double, double>> partPos```` lists the positions of each parton. To add a parton, the user adds a new pair of x and y values to ````partPos````. ````partNum```` is updated by calculating the new size of the vector.
