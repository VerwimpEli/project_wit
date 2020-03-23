#include "BoundaryCondition.cpp"
#include "FVM.cpp"
#include "adjoint.cpp"
#include  <Eigen/Sparse>
#include <vector>

/**
 * Temporary version of Matlab scale function applied to -1 vector. VW might have to be VH, not sure.
 */
void build_rhs(Eigen::VectorXd & v, int VW){
    double value = -1;
    for (int i = 0; i < v.rows(); i++) {
        if (i < VW){ value *= 0.5; }
        else if (i > v.rows() - VW - 1){ value *= 0.5; }

        if (i % VW == 0){ value *= 0.5; }
        else if ((i + 1) % VW == 0) {value *= 0.5; }

        v(i) = value;
        value = -1;
    }
}

int main(){

    double H = 1.0;
    double W = 1.0;
    int VW = 25;
    int VH = 25;
    double Q = 200.0;
    double Cmet = 65.0;
    double Cpla = 0.2;
    double M = 0.4;
    int p = 5;

    std::vector<double> mat(VW*VH, 0.2);


    // Boundary Conditions

    BoundarySegment b1 = BoundarySegment(NEUMANN, 0, 1, 0);
    BoundaryCondition BcBottom(b1);

    BoundarySegment r1 = BoundarySegment(NEUMANN, 0, 0.3, 0);
    BoundarySegment r2 = BoundarySegment(DIRICHLET, 0.3, 0.7, 293);
    BoundarySegment r3 = BoundarySegment(NEUMANN, 0.7, 1, 0);

    std::vector<BoundarySegment> SegVecRight({r1, r2, r3});
    BoundaryCondition BcRight(SegVecRight);

    BoundaryCondition BcTop(b1);
    BoundaryCondition BcLeft(SegVecRight);

    // Solve

    std::vector<double> Sol(VW*VH);
    Eigen::SparseMatrix<double> K(VW * VH, VW* VH);

    FVM<double> fvm(H, W, VW, VH, Q, Cmet, Cpla, p, BcBottom, BcRight, BcTop, BcLeft);
    AdjointGradient ag(H, W, VW, VH, M, Q, Cmet, Cpla, p);

    fvm(mat, Sol, K);

    K = K.transpose();
    Eigen::VectorXd rhs(VW * VH);
    build_rhs(rhs, VW);

    Eigen::VectorXd L(VW * VH);
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::NaturalOrdering<int>> solver;
    solver.compute(K);
    L = solver.solve(rhs);

    std::vector<double> AG(VW * VH);
    ag(mat, L, Sol, AG);

    std::cout << "\n\n---------- Solution ----------\n\n";
    Print(Sol);
    std::cout << "\n\n------------ AG --------------\n\n";
    Print(AG);

    return 0;
}