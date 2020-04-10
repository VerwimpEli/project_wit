CXX = g++
CXXFLAGS = -O3 -fopenmp -g
INCLUDE = -I./dependencies/eigen/ 
INCLUDESS = -I./dependencies/Linux/SuiteSparse/include/
LIBINCLUDE = -L ./dependencies/Linux/SuiteSparse/lib
LDFLAGS = -lumfpack -lamd -lsuitesparseconfig -lm

main: main.o MMASolver.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o apps/main main.o MMASolver.o

main.o: src/main.cpp 
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c src/main.cpp

MMASolver.o: src/MMASolver.cpp
	$(CXX) $(CXXFLAGS) -c src/MMASolver.cpp

main_umf: main.o MMASolver.o
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(INCLUDESS) $(LIBINCLUDE)  -o apps/main_umf main.o MMASolver.o $(LDFLAGS)  

clean:
	-rm *.o

# /**
#  *  Suitesparse solver, fastest option. Requires: Lapack, Openblas, SuiteSparse and OpenMP
#  *  g++ -O3 -fopenmp -I ./eigen/ -I ./SuiteSparse-master/include/ main.cpp MMASolver.cpp -lumfpack -o main.o
#  *
#  *  Normal solver, only requires EIGEN (but has speedup with lapack and openblas)
#  *  g++ -O3 -I ./eigen/ main.cpp MMASolver.cpp -o main.o
#  */
