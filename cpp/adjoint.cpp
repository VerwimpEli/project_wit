#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>

class AdjointGradient {
    public: 
        std::vector<double> AdjointGradient::AdjointGradient(int VW, int VH, float M, float Q, float Cmet,
        float Cpla, int p, BoundaryCondition BC[])
            : VW_(VW)
            , VH_(VH)
            , M_(M)
            , Q_(Q)
            , Cmet_(Cmet)
            , Cpla_(Cpla)
            , p_(p)
            , BC_(BC)
            ;
}