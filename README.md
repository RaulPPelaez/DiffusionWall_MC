# DiffusionWall_MC  

Raul P. Pelaez 2018.   
  
  An UAMMD[1] Monte Carlo simulation of Debye type blobs confined near a wall with a Debye repulsion.  
  It uses the MonteCarloNVT_Anderson module by Pablo Ibañez[2].  
  
  ### USAGE  
  
  ## HOW TO COMPILE DiffusionWall_MC.cu  
  
  Compile with nvcc, you need to tell it where UAMMD and HydroGrid[3] are (HydroGrid has to be compiled already) and the flag --expt-relaxed-constexpr is necessary.  
  You will need to link HydroGrid in some way. I usually pass -lCallHydroGrid and then copy libCallHydroGrid.so to the same folder as the executable, it ight be different for your system.  
  Refer to [5] for more compilation options and details. You should have recieved a Makefile along this file though.  
  
  
  ## HOW TO RUN A SIMULATION  
  
  This code needs the following files to be in the same folder in order to run:  
  
  -data.main            -> A file containing the simulation parameters, you should have received an example data.main along this file.  
  -hydroGridOptions.nml -> A HydroGrid configuration file. See [3] for more info.  
  
  After compiling, just execute the binary. You can pass the --device X flag to select a certain GPU.  
  
  
  A simulation is already prepared under simulations/DynamicStructureFactor
  
  Just go into this folder and run bash run.bash. This script will take care of compiling and running the simulation.
  DynamicStructureFactor is configured to run a simulation, compute dynamic structure factors with HG, radial distribution function and probability distribution in z.
  After it is done, you will find a DynamicStructureFactor/tools/sim.results folder with all results.
  
  ### SYSTEM DETAILS   
  
  Particles interact via a repulsive Debye potential.  
  
  Particles interact with two walls in the xy plane.    
  The walls are:  
  1- At Lz/2 (infinite wall)  
  2- At -Lz/2(infinite wall + Debye repulsion)  
  
  
  The system is technically unbounded in the z direction. The two inifinite zwalls avoid periodic artifacts like particles seeing each other through the floor/ceiling of the box. The box is increased by blob_cutOff in the z direction to avoid particles near one of the walls interact with the other due to PBC.  
  
  ### PARAMETERS AND INITIAL CONFIGURATION  
  
  All parameters are read from a data.main file. You should have recieved an example reference one with this file.  
  You can refer to the function readParameters to see a list of available parameters.  
  
  The initial configuration of particles will be random if a coordinates file is not provided in data.main  
  The format of the coordinates file provided in the data.main is:  
  ----  
  Nparticles  
  x y z  
  ...  
  ----  
  A fourth column might be included in the file that will be interpreted as species by HydroGrid. If only three columns are available all particle are assumed to have species 0.  
  
  ### HYDROGRID ANALYSIS  
  
  This code is able to call HydroGrid[3] to process the simulation results on the fly.  
  
  HydroGrid is called every samplefreq to update itself according to hydroGridOptions.nml (which must be present at the current folder) and write to disk. You should have received an example .nml with this file, see [3] for more information about HydroGrid and its input.  
  
  HydroGrid is called through the UAMMD wrapper class HydroGrid in utils/HydroGrid.cuh, you can see more details there.  
  
  ### OUTPUT  
  
  Information about the current simulation will be output at runtime through stdout and stderr,  
  including all parameters used and some useful information.   
  
  Besides HydroGrid output, this code will output the particle positions every savefreq with the following format:  
  ----  
  #time1  
  x y z radius color  
  ...  
  #time2  
  ...  
  ----  
  The positions will not be folded to the primary simulation box.  
  
  You can visualize the particle results with superpunto[4]. HG will output whatever was asked in the .nml.  
  
  
  ### References:   
  [1] https://github.com/RaulPPelaez/UAMMD/wiki  
  [2] https://github.com/PabloIbannez/UAMMD/blob/master/src/Integrator/MonteCarlo/MonteCarloNVT_Anderson.cuh  
  [3] https://github.com/stochasticHydroTools/HydroGrid  
  [4] https://github.com/RaulPPelaez/superpunto  
  [5] https://github.com/RaulPPelaez/UAMMD/wiki/Compiling-UAMMD  
