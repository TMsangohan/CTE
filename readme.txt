In this folder the CTE code (status 11/02/2015) including
- betastar levelling
- stochastic cooling 
can be found. 

The files 
collider_time_evolution_cool.f90
collider_time_evoluiton_common.f90
contain the program code.

collider_time_evolution_cool.in 
is the required input file.

The madx twiss table 
LHCB1-CollIons2011-6.503.tfs
was put here as an example optics to show the structure of the required file.

To compile the code run the following command
gfortran -fdefault-real-8 -fdefault-double-8 -fno-align-commons -o collider_time_evolution collider_time_evolution.f90 

the executable
collider_time_evolution
is produced and is to be called as
.\collider_time_evolution 
to run the simulation with all the above files in the same folder. The output will also be written to the folder from that the code was run.

All .out files were produced in a test run of the code.

The folder ExampleBatchJobSub contains an example how to submitt simulations are batch jobs. For that run
.\translateDos2Unix
This file will translate the batchit and .in files produced under Windows for example in Mathematica to the unix line endings.
Makes the runit file an executable and calls it. runit contains all jobs to be submitted in one go.
Each job requires an .in and one batchit file. The batchit file contains information of where to find the simulation code and input files, lattice files etc., it starts the simulaiton and copies the output back to a choosen location.


Detailed information about the structure and varables of the input file, can be found in the appendix of my PhD thesis:
Michaela Schaumann
Heavy-Ion Perfomance of the LHC and Future Colliders.
