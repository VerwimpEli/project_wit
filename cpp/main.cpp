#include "MMASolver.h"
#include "BoundaryCondition.cpp"
#include "FVM.cpp"
#include "adjoint.cpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <numeric>

template<typename T>
void scale(T & v, int VW, int size){
    for (int i = 0; i < size; i++) {
        if (i < VW){ v[i] *= 0.5; }
        else if (i > size - VW - 1){ v[i] *= 0.5; }

        if (i % VW == 0){ v[i] *= 0.5; }
        else if ((i + 1) % VW == 0) {v[i] *= 0.5; }
    }
}

class HeatEq
{
public:
    HeatEq(double H, double W, int VW, int VH, double Q, double Cmet, double Cpla, int p, double M,
           BoundaryCondition BC0, BoundaryCondition BC1, BoundaryCondition BC2, BoundaryCondition BC3)
            : H_(H), W_(W), Cmet_(Cmet), Cpla_(Cpla), VW_(VW), VH_(VH), p_(p), Q_(Q), M_(M),
              BC0_(BC0), BC1_(BC1), BC2_(BC2), BC3_(BC3),
              fvm(H, W, VW, VH, Q, Cmet, Cpla, p, BC0, BC1, BC2, BC3),
              ag(W, H, VW, VH, M, Q, Cmet, Cpla, p)
    {}

    void operator() (std::vector<double> & v, double & f, std::vector<double> & df,
                     double & g, std::vector<double> & dg){

        std::vector<double> sol(VW_ * VH_);
        Eigen::SparseMatrix<double> K(VW_ * VH_, VW_ * VH_);

        fvm(v, sol, K);
        Print(sol);
        f = std::accumulate(sol.begin(), sol.end(), 0.0);

        K = K.transpose();
        Eigen::VectorXd rhs(VW_ * VH_);
        for (int i = 0; i < VW_ * VH_; i++){
            rhs(i) = -1;
        }
        scale(rhs, VW_, rhs.rows());

        Eigen::VectorXd L(VW_ * VH_);
        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::NaturalOrdering<int>> solver;
        solver.compute(K);
        L = solver.solve(rhs);

        ag(v, L, sol, df);

        std::vector<double> vCopy(v);
        scale(vCopy, VW_, vCopy.size());
        g = std::accumulate(vCopy.begin(), vCopy.end(), 0.0) - M_ * VW_ * VH_;

        std::fill(dg.begin(), dg.end(), 1.0);
        scale(dg, VW_, dg.size());
    }

    void solve_heat_eq(std::vector<double> & v, std::vector<double> & t){
        Eigen::SparseMatrix<double> K(VW_ * VH_, VW_ * VH_);
        fvm(v, t, K);
    }



private:
    double H_, W_, Cmet_, Cpla_, Q_, M_;
    int VW_, VH_, p_;
    BoundaryCondition BC0_, BC1_, BC2_, BC3_;
    FVM<double> fvm;
    AdjointGradient ag;

};

template <typename T>
T Min(T d1, T d2) {
    return d1<d2 ? d1 : d2;
}

template <typename T>
T Max(T d1, T d2) {
    return d1>d2 ? d1 : d2;
}

double Abs(double d1) {
    return d1>0 ? d1 : -1.0*d1;
}


int main(int argc, char *argv[]) {

    double H = 1.0;
    double W = 1.0;
    int VW = 10;
    int VH = 10;
    double Q = 200.0;
    double Cmet = 65.0;
    double Cpla = 0.2;
    double M = 0.4;
    int p = 1;

    std::vector<double> v(VW*VH, 0.1);

    BoundarySegment b1 = BoundarySegment(NEUMANN, 0, 1, 0);
    BoundaryCondition BcBottom(b1);

    BoundarySegment r1 = BoundarySegment(NEUMANN, 0, 0.3, 0);
    BoundarySegment r2 = BoundarySegment(DIRICHLET, 0.3, 0.7, 293);
    BoundarySegment r3 = BoundarySegment(NEUMANN, 0.7, 1, 0);

    std::vector<BoundarySegment> SegVecRight({r1, r2, r3});
    BoundaryCondition BcRight(SegVecRight);

    BoundaryCondition BcTop(b1);
    BoundaryCondition BcLeft(SegVecRight);

	int n = VW * VH;
	int m = 1;
	double movlim = 1.0;
	int max_iter = 100;

    std::vector<double> v_old(v);
    double f;
    std::vector<double> df(n);
    double g;
    std::vector<double> dg(n);
    std::vector<double> v_min(n, 0);
    std::vector<double> v_max(n, 1);

	MMASolver *mma = new MMASolver(n,m);
	HeatEq heateq(H, W, VW, VH, Q, Cmet, Cpla, p, M, BcBottom, BcRight, BcTop, BcLeft);

	double ch = 1.0;
	int itr = 0;
	while (itr < max_iter) {
		itr++;

		heateq(v, f, df, g, dg);

//        std::cout << "\n\n---------- f ----------\n\n";
//        std::cout << f << std::endl;
//        std::cout << "\n\n---------- g ----------\n\n";
//        std::cout << g << std::endl;
//        std::cout << "\n\n---------- df ----------\n\n";
//		Print(df);
//		std::cout << "\n\n--------- dg --------\n\n";
//		Print(df);

		// Set outer move limits
		for (int i=0;i<n;i++) {
			v_max[i] = Min(v_max[i], v[i] + movlim);
			v_min[i] = Max(v_min[i], v[i] - movlim);
		}

		// Call the update method
		double* v_arr = &v[0]; // Optimization works with arrays. Might change later.
        double* df_arr = &df[0];
        double* g_arr = &g;
        double* dg_arr = &dg[0];
        double* v_min_arr = &v_min[0];
        double* v_max_arr = &v_max[0];

        mma->Update(v_arr,df_arr,g_arr,dg_arr,v_min_arr,v_max_arr);

		// Compute infnorm on design change
		ch = 0.0;
		for (int i=0;i<n;i++) {
			ch = Max(ch,Abs(v[i]-v_old[i]));
			v_old[i] = v[i];
		}

		// Print to screen
		printf("it.: %d, obj.: %f, g.: %f, ch.: %f \n",itr, f, g, ch);
	}

	std::cout << "v:";
	for (int i=0;i<n;i++) {
		std::cout << " " << v[i];
	}
	std::cout << std::endl;

	std::vector<double> t(n);
	heateq.solve_heat_eq(v, t);

	std::ofstream file;
	file.open("test.out");
	for (int i = 0; i < n -1; i++) {
	    file << v[i] << ",";
	}
	file << v[n - 1] << "\n";

    for (int i = 0; i < n -1; i++) {
        file << t[i] << ",";
    }
    file << t[n - 1] << "\n";

	file.close();

	// Deallocate
	delete mma;
	return 0;

}

