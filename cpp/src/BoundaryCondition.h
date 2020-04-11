#ifndef CPP_BOUNDARYCONDITION_H
#define CPP_BOUNDARYCONDITION_H

#define  DIRICHLET 0
#define  NEUMANN 1

#include <vector>
class BoundarySegment {
public:
    BoundarySegment(int type, float start, float stop, float value);
    BoundarySegment(BoundarySegment const & bs);
    BoundarySegment();

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

class BoundaryCondition {
public:
    BoundaryCondition(std::vector<BoundarySegment> const s);
    BoundaryCondition(BoundarySegment const s);
    std::vector<BoundarySegment> GetSegments() const { return segments_; }
    BoundarySegment GetStart() const { return start; }
    BoundarySegment GetStop() const { return stop; }
private:
    BoundarySegment start;
    BoundarySegment stop;
    std::vector<BoundarySegment> segments_;
};
#endif //CPP_BOUNDARYCONDITION_H
