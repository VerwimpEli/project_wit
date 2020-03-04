//
// Created by eli on 04/03/2020.
//
#include <assert.h>
#include <iostream>

#ifndef DIRICHLET
#define  DIRICHLET 0
#endif

#ifndef NEUMANN
#define  NEUMANN 1
#endif

class BoundaryCondition
{
public:
    BoundaryCondition(int type, int start, int stop, int value)
    : type_(type), start_(start), stop_(stop), value_(value) {
        assert(type == 0 || type == 1);
    }

    int type() { return type_; }
    int start() { return start_; }
    int stop() { return stop_; }
    int value() { return value_; }

private:
    int type_;
    int start_;
    int stop_;
    int value_;
};

/**
 * Test/Example of BoundaryCondition
 */
int main(){

    BoundaryCondition b = BoundaryCondition(DIRICHLET, 0.4, 0.6, 293);
    std::cout << b.value() << std::endl;
    std::cout << b.type() << std::endl;

    return 0;
}