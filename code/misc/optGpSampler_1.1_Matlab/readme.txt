==================== optGpSampler version 1.1 ==================
Installing and testing optGpSampler for MATLAB.

optGpSampler has pre-compiled connectors for MATLAB (Win/Linux). It requires any of the following linear programming (LP) solvers:
a) IBM ILOG CPLEX 12.6 or later
b) GUROBI 5.6 or later
c) GLPK 4.53 or later

=== Running optGpSampler on Linux ===
1) Add the "linux_lib" subdirectory to your LD_LIBRARY_PATH
2) Make sure that the path to the dynamic library of your solver is in your LD_LIBRARY_PATH.
- For CPLEX: the path to "cplex1260.so"
- For GUROBI: the path to "gurobi56.so"
- For GLPK: glpk.so (shipped within the linux_lib subdirectory of the sampler)
3) Do not rename or remove the directory "linux_lib". Libraries are loaded dynamically at runtime and are expected to be in this location.

=== Running optGpSampler on Windows ===
1) Make sure the path to your solver is in your system path variable (this is usually the default if you install the solver).
2) Load a model and run optGpSampler.m

=== Workaround for other versions ===
If you do not have the exact version of the LP library, it usually works when you copy the library
e.g. gurobi55.dll and rename it to the correct version. (in this case gurobi56.dll). 
For Windows, copy the renamed file to the sampler folder
For Linux, make sure the renamed file is in the LD_LIBRARY_PATH


== CLEAN-UP (OPTIONAL) == 
The sampler comes with files for both Linux and Windows, and connectors to all supported solvers. If you like to clean up files you don't use, you can do the following:
LINUX: remove all *.dll files, all the *.so files for solvers you don't use and the *.mexW64 files (mind that you don't remove the *.mexA64 files)
WINDOWS: remove the linux_lib directory, all *.dll files to solvers you don't use and the *.mexA64 files (mind that you don't remove the *.mexW64 files)