# Debris Cloud Dynamics Rebound Package

This package allows for the simulation of debris clouds using N-body integrator Rebound. For specifics on the construction of this package, please reference Larson and Sarid (2021). This document should act as a manual for using this package. 

Necessary files/modules necessary to use package: 

-Python version of Rebound

-Python (any version works)

-Main script (binarytest.py)

-Input files (bodyparams.py; integratorparams.py)

-Effects files (additionaleffects.py; basicsolarsystem1015.py)

-Plotting script (plotejecta2.py)
 
	
	
# Simulation Set-up

All files must be saved in the same directory. The imput files bodyparams.py and bodyparams.py adjust the initial parameters of the simulation. The file bodyparams-sample.txt is a sample of the bodyparams.py input file with descriptions of each input variable. Similarly, integratorparams-sample.txt describes the input variables for integratorparams.py.

After the input files are set up with the initial parameters of the simulation, the simulation runs by calling the binarytest.py script with python in the terminal like so:

username:~/project_directory$ python binarytest.py

The output files with particle data at each time step will be saved to the project directory. It is possible to write a new plotting function to present the data; however, we have also included plotejecta2.py as a very simplistic plotting function that plots the particles in the x-y plane and the x-z plane. These two views may not be ideal for the scope of the user's study in which case we encourage them to use plotejecta2.py as a baseline from which to construct a new plotting function or to simply ignore plotejecta2.py and proceed with their own specific data analysis. For ease of use, plotejecta2-sample.txt explains the set-up of the script in order to plot the data. To run the plotejecta2.py script, call the file with python in the terminal:

username:~/project_directory$ python plotejecta2.py

We recommend plotting the data files in the same simulation directory in which the simulation was run in order to keep track of which plots are associated with which simulation. Later, after the plots are made, if the user wishes to save the plots in a separate directory, they can easily copy the plots to a separate directory (preferably named to reflect to which simulation those plots pertain).


# Function Descriptions: basicsolarsystem1015.py



# Function Descriptions: additionaleffects.py
