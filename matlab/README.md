## Project WIT: Matlab

## Main script
The `main.m` script will solve the main assingment. Most significant options that 
maybe want to be changed are `VW` and `VH` which denote the gridsize in the width and heigth.
Also `maxiter` can be adapted to change the number of optimization iterations.

The main script will plot two figures. Figure 1 will have the final solution and figure 2 will
have the inital temperature and the optimized temperature.
 
## Tests

Tests are comprised of 2 parts. The first part is matlab code, with the purpose to prove the correctness of the matlab implementation. The second part is part matlab, part c++ and has the goal to prove that both give the same result. This would mean that the c++ code also delivers a correct result.

FIRST PART : The tests are split in to files. The first test the correctness of the Finite volume method. These tests are in the 'Harmonic_TestCase_FVM'. More explaintion of these tests is included in the comments of this file. The second test is 'Harmonic_TestCase_adjoint' and compares the gradient computed with the adjoint method with a finite difference approximation. 

SECOND PART : in matlab we have 'Harmonic_TestCase_adjoint' and in c++ we use 'fvm_adjoint_test.cpp'. Both generate the temperature solution and the gradient associated with that solution. In order to avoid the requirement of matlab, a reference solution is provided.   

