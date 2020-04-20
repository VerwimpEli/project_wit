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
The Makefile in this folder is able to build all programs (on Linux). The default compiler is g++ and
the default options include -O3 and -Wall. These can be changed through the variables CXX and CXXFLAGS.
 There are a couple options:

- make: default make will build the main script
- make test: Build & compile the test files
- make time_main: Build & compile the timed main script
- make main_umf: Build & compile the main script, with SuiteSparse's UmfLUSolver. Suitesparse is included, but 
OpenBlas and Lapack libraries should be installed on the system. SuiteSparse is dynamically linked so when running
the system should be able to find it. (i.e. add to LD_LIBRARY_PATH in Linux). See below for exhausive instructions.
- make time_main_umf: Similar to time_main but with UmfLUSolver. See above for comments on how to use the UmfLUSolver.
- make all: make all of the above
- make clean: delete object files


Step by step linux (ubuntu) installation to make main_umf:

_These instructions also work on the linux for windows subsystem._
- Install openblas: `sudo apt install libopenblas-dev`
- Install lapack: `sudo apt install liblapack-dev`
- Add SuiteSparse to library path: 
`export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:INSTALL_DIR/project_wit/cpp/dependencies/Linux/SuiteSparse/lib"`.
 Where INSTALL_DIR is replaced by your installation path.
- Make script: `make main_umf`

#### Compile & build on windows
_Unless somebody can compile SuiteSparse on windows, the fast UmfLUSolver is linux only. Test fails can be compiled in the same way as the main program.
Replace CXX by favorite compiler._
- `CXX -O3 -Wall -I./dependencies/eigen/ -I./src/ -o apps/main_win src/main.cpp src/MMASolver.cpp src/util.cpp src/BoundaryCondition.cpp` 

### Main programs
_Note: run programs from within apps folder. Else the results won't be stored correctly_

The **main** program will solve the heat topology problem with the given configuration. Two variables can be passed through the
command line: `max_iterations` and `size`. These are respectively the maximum number of iterations the MMA Solver will 
perform and the number of design variables in the width and height of the domain. So there are size² design variables. _(ex: ./main 500 256, for
500 iterations on 256 by 256 domain_). Defaults to 500 iterations and 32x32 domain. Other options can be changed in `main.cpp`.


The result will be stored in ./results/ in an .out file, with a datetime stamp in its name. This file contains the last
v vector as well as the temperature vector. These can be plot with the python script in this folder: `python plot.py -f [file] -u [size]`

The **main_umf** is the same as the main program, but with SuiteSparse's fast UmfLUSolver. 
On smaller grids there is an overhead, but for grids larger than 256 this solver should
be used. Details are the same as the main program.

The **time_main_(umf)** programs run one iteration of the main programs and record the time
needed to solve one FVM calculation, one adjoint calculation and one MMA optimization step. 
Have the same options as main programs, but the first (`max_iterations`) is ignored.

### Tests

Tests are comprised of 2 parts. The first part is matlab code, with the purpose to prove the correctness of the matlab implementation. The second part is part matlab, part c++ and has the goal to prove that both give the same result. This would mean that the c++ code also delivers a correct result.

FIRST PART : The tests are split in to files. The first test the correctness of the Finite volume method.
These tests are in the 'TestCase_FVM'. More explaintion of these tests is included in the comments of this file.
The second test is 'TestCase_adjoint' and compares the gradient computed with the adjoint method with a finite difference approximation. 

SECOND PART : in matlab we have 'TestCase_adjoint' and in c++ we use 'fvm_adjoint_test.cpp'.
Both generate the temperature solution and the gradient associated with that solution.
In order to avoid the requirement of matlab, a reference solution is provided. 
