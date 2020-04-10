CXX = g++
CXXFLAGS = -O3 -Wall
UMFFLAG = -DUmfLUSolver
INCLUDE = -I./dependencies/eigen/ -I./src/
INCLUDESS = -I./dependencies/Linux/SuiteSparse/include/
LIBINCLUDE = -L ./dependencies/Linux/SuiteSparse/lib
LDFLAGS = -lumfpack -lamd -lsuitesparseconfig -lm

main: main.o MMASolver.o util.o BoundaryCondition.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o apps/main main.o MMASolver.o util.o BoundaryCondition.o 

main.o: src/main.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c src/main.cpp

MMASolver.o: src/MMASolver.cpp
	$(CXX) $(CXXFLAGS) -c src/MMASolver.cpp

util.o: src/util.cpp
	$(CXX) $(CXXFLAGS) -c src/util.cpp

BoundaryCondition.o: src/BoundaryCondition.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c src/BoundaryCondition.cpp

main_umf: main_umf.o MMASolver.o util.o BoundaryCondition.o
	$(CXX) $(CXXFLAGS) $(UMFFLAG) $(INCLUDE) $(INCLUDESS) $(LIBINCLUDE)  -o apps/main_umf main_umf.o MMASolver.o util.o BoundaryCondition.o $(LDFLAGS)  

main_umf.o: src/main.cpp 
	$(CXX) $(CXXFLAGS) $(UMFFLAG) $(INCLUDE) $(INCLUDESS) -c src/main.cpp -o main_umf.o 

fvm_adjoint_test: util.o BoundaryCondition.o 
	$(CXX) $(CXXFLAGS) $(INCLUDE) test/fvm_adjoint_test.cpp -o apps/fvm_adjoint_test util.o BoundaryCondition.o

time_main: time_main.o MMASolver.o util.o BoundaryCondition.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) -DTIME -o apps/time_main time_main.o MMASolver.o util.o BoundaryCondition.o 

time_main.o: src/main.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -DTIME -c src/main.cpp -o time_main.o

time_main_umf.o: src/main.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(UMFFLAG) $(INCLUDESS) -DTIME -c src/main.cpp -o time_main_umf.o

time_main_umf: time_main_umf.o MMASolver.o util.o BoundaryCondition.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(UMFFLAG) $(INCLUDESS) $(LIBINCLUDE) -DTIME -o apps/time_main_umf time_main_umf.o MMASolver.o util.o BoundaryCondition.o $(LDFLAGS)

test: fvm_adjoint_test


all: main_umf main fvm_adjoint_test

clean:
	-rm *.o

# /**
#  *  Suitesparse solver, fastest option. Requires: Lapack, Openblas, SuiteSparse and OpenMP
#  *  g++ -O3 -fopenmp -I ./eigen/ -I ./SuiteSparse-master/include/ main.cpp MMASolver.cpp -lumfpack -o main.o
#  *
#  *  Normal solver, only requires EIGEN (but has speedup with lapack and openblas)
#  *  g++ -O3 -I ./eigen/ main.cpp MMASolver.cpp -o main.o
#  */
