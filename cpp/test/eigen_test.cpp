#include <iostream>
#include <Eigen/Dense>

/**
 * CD into cpp folder
 * GET EIGEN: git clone https://gitlab.com/libeigen/eigen.git
 * COMPILE: g++ -I ./eigen/ eigen_test.cpp -o eigen_test.o
 */

using Eigen::MatrixXd;
int main()
{
    MatrixXd m(2,2);
    m(0,0) = 3;
    m(1,0) = 2.5;
    m(0,1) = -1;
    m(1,1) = m(1,0) + m(0,1);
    std::cout << m << std::endl;
}