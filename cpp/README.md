## CPP: project WIT
### Project structure
- ./apps: compiled and build programs will be stored here.
- ./dependencies: contains (some) of the dependencies for the programs, see
compile and build for more information.
- ./results: Output of main script will appear here. There's also a 
python script here to plot the results
- ./src: Containts the source files of the project. This is the cpp code 
and its header files.
- ./test: Source files for the tests and test scripts for cpp are stored here.
### Compile & build
The Makefile in this folder are able to build all script (currently only linux). The default compiler is g++ and
the default options include -O3 and -Wall. These can be changed through the variables CXX and CXXFLAGS.
 There are a couple options:

- make: default make will build the main script
- make test: Build & compile the test files
- make time_main: Build & compile the timed main script
- make main_umf: Build & compile the main script, with SuiteSparse's UmfLUSolver. Suitesparse is included, but 
OpenBlas and Lapack libraries should be installed on the system. SuiteSparse is dynamically linked so when running
the system should be able to find it. (i.e. add to LD_LIBRARY_PATH in Linux).
- make time_main_umf: Similar to time_main but with UmfLUSolver. See above for comments on how to use the UmfLUSolver.
- make clean: delete object files
### Tests

### Main