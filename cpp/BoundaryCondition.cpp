//
// Created by eli on 04/03/2020.
//
#include <assert.h>
#include <iostream>
#include <vector>

#ifndef DIRICHLET
#define  DIRICHLET 0
#endif

#ifndef NEUMANN
#define  NEUMANN 1
#endif

class BoundarySegment
{
public:
    BoundarySegment(int type, float start, float stop, float value)
    : type_(type), start_(start), stop_(stop), value_(value) {
        assert(type == 0 || type == 1);
    }

    BoundarySegment(BoundarySegment const & bs)
    : type_(bs.type_), start_(bs.start_), stop_(bs.stop_), value_(bs.value_) {}

    BoundarySegment()
    : type_(DIRICHLET), start_(0), stop_(1), value_(273) {}

    int type() const { return type_; }
    float start() const { return start_; }
    float stop() const { return stop_; }
    float value() const { return value_; }

private:
    int type_;
    float start_;
    float stop_;
    float value_;
};


class BoundaryCondition
{
public:
    BoundaryCondition(std::vector<BoundarySegment> const s)
    : segments_(s) {

        for (BoundarySegment seg: s){
            if (seg.start() == 0){
                start = BoundarySegment(seg.type(), 0, 0, seg.value());
            }

            if (seg.stop() == 1){
                stop = BoundarySegment(seg.type(), 1, 1, seg.value());
            }
        }
    }

    BoundaryCondition(BoundarySegment const s)
    : segments_({s}) {
        assert(s.start() == 0);
        assert(s.stop() == 1);
        start = BoundarySegment(s.type(), 0, 0, s.value());
        stop = BoundarySegment(s.type(), 1, 1, s.value());
    }

    std::vector<BoundarySegment> GetSegments() { return segments_; }
    BoundarySegment GetStart(){ return start; }
    BoundarySegment GetStop(){ return stop; }

private:
    BoundarySegment start;
    BoundarySegment stop;
    std::vector<BoundarySegment> segments_;
};


/**
 * Test/Example of BoundaryCondition
 * Don't outcomment this if using in other files.
 */
//int main(){
//
//    BoundarySegment a = BoundarySegment(NEUMANN, 0.6, 1, 0);
//    BoundarySegment b = BoundarySegment(DIRICHLET, 0.4, 0.6, 293);
//    BoundarySegment c = BoundarySegment(NEUMANN, 0, 0.4, 5);
//    std::vector<BoundarySegment> SegVec({a, b, c});
//    BoundaryCondition bc(SegVec);
//
//    for (BoundarySegment seg : bc.GetSegments()) {
//        std::cout << seg.type() << std::endl;
//    }
//
//    std::cout << "Start: " << bc.GetStart().type() << std::endl;
//    std::cout << "Start: " << bc.GetStart().start() << std::endl;
//    std::cout << "Start: " << bc.GetStart().stop() << std::endl;
//    std::cout << "Start: " << bc.GetStart().value() << std::endl;
//
//    std::cout << "Stop: " << bc.GetStop().type() << std::endl;
//    std::cout << "Stop: " << bc.GetStop().start() << std::endl;
//    std::cout << "Stop: " << bc.GetStop().stop() << std::endl;
//    std::cout << "Stop: " << bc.GetStop().value() << std::endl;
//    return 0;
//}