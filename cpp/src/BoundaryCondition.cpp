#include <assert.h>
#include "BoundaryCondition.h"
#include <iostream>
#include <vector>

#ifndef DIRICHLET
#define  DIRICHLET 0
#endif

#ifndef NEUMANN
#define  NEUMANN 1
#endif


BoundarySegment::BoundarySegment(int type, float start, float stop, float value)
    : type_(type), start_(start), stop_(stop), value_(value) {
          assert(type == 0 || type == 1);
    }

BoundarySegment::BoundarySegment(BoundarySegment const & bs)
    : type_(bs.type_), start_(bs.start_), stop_(bs.stop_), value_(bs.value_) {}

BoundarySegment::BoundarySegment()
    : type_(DIRICHLET), start_(0), stop_(1), value_(273) {}

BoundaryCondition::BoundaryCondition(std::vector<BoundarySegment> const s)
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

BoundaryCondition::BoundaryCondition(BoundarySegment const s)
            : segments_({s}) {
    assert(s.start() == 0);
    assert(s.stop() == 1);
    start = BoundarySegment(s.type(), 0, 0, s.value());
    stop = BoundarySegment(s.type(), 1, 1, s.value());
}





//class BoundarySegment
//{
//public:
//    BoundarySegment(int type, float start, float stop, float value)
//    : type_(type), start_(start), stop_(stop), value_(value) {
//        assert(type == 0 || type == 1);
//    }
//
//    BoundarySegment(BoundarySegment const & bs)
//    : type_(bs.type_), start_(bs.start_), stop_(bs.stop_), value_(bs.value_) {}
//
//    BoundarySegment()
//    : type_(DIRICHLET), start_(0), stop_(1), value_(273) {}
//
//    int type() const { return type_; }
//    float start() const { return start_; }
//    float stop() const { return stop_; }
//    float value() const { return value_; }
//
//private:
//    int type_;
//    float start_;
//    float stop_;
//    float value_;
//};
//
//
//class BoundaryCondition
//{
//public:
//    BoundaryCondition(std::vector<BoundarySegment> const s)
//    : segments_(s) {
//
//        for (BoundarySegment seg: s){
//            if (seg.start() == 0){
//                start = BoundarySegment(seg.type(), 0, 0, seg.value());
//            }
//
//            if (seg.stop() == 1){
//                stop = BoundarySegment(seg.type(), 1, 1, seg.value());
//            }
//        }
//    }
//
//    BoundaryCondition(BoundarySegment const s)
//    : segments_({s}) {
//        assert(s.start() == 0);
//        assert(s.stop() == 1);
//        start = BoundarySegment(s.type(), 0, 0, s.value());
//        stop = BoundarySegment(s.type(), 1, 1, s.value());
//    }
//
//    std::vector<BoundarySegment> GetSegments() const { return segments_; }
//    BoundarySegment GetStart() const { return start; }
//    BoundarySegment GetStop() const { return stop; }
//
//private:
//    BoundarySegment start;
//    BoundarySegment stop;
//    std::vector<BoundarySegment> segments_;
//};