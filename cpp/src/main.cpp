#define TIME true           // Set to false to turn of timings
#define UmfLUSolver false    // Set to false for simple solver. (VW <= 256)
#define MT true             // Multithread on

#include <chrono>
#include "MMASolver.h"
#include <iostream>
#include "util.cpp"
#include "BoundaryCondition.cpp"
#if UmfLUSolver
#include <Eigen/UmfPackSupport>
#endif
#include "FVM.cpp"
#include "adjoint.cpp"
#include "MMASolver.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <numeric>
#include <ctime>
#if MT
#include <omp.h>
#endif


/**
 *  Suitesparse solver, fastest option. Requires: Lapack, Openblas, SuiteSparse and OpenMP
 *  g++ -O3 -fopenmp -I ./eigen/ -I ./SuiteSparse-master/include/ main.cpp MMASolver.cpp -lumfpack -o main.o
 *
 *  Normal solver, only requires EIGEN (but has speedup with lapack and openblas)
 *  g++ -O3 -I ./eigen/ main.cpp MMASolver.cpp -o main.o
 */

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
        Eigen::VectorXd L(VW_ * VH_);

        auto t_start = std::chrono::system_clock::now();

        fvm(v, sol, K, L);

        #if TIME
            auto t_end = std::chrono::system_clock::now();
            std::chrono::duration<double> diff = t_end-t_start;
            std::cout << "FVM took: " << diff.count() << " s" << std::endl;
        #endif

        t_start = std::chrono::system_clock::now();

        ag(v, L, sol, df);

        #if TIME
            t_end = std::chrono::system_clock::now();
            diff = t_end-t_start;
            std::cout << "ADJ took: " << diff.count() << " s" << std::endl;
        #endif

        f = std::accumulate(sol.begin(), sol.end(), 0.0);

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

    void update_p(int p){
        std::cout << "Updated p to: " << p << std::endl;
        p_ = p;
        fvm.update_p(p);
        ag.update_p(p);
    }


private:
    double H_, W_, Cmet_, Cpla_, Q_, M_;
    int VW_, VH_, p_;
    BoundaryCondition BC0_, BC1_, BC2_, BC3_;
    FVM<double> fvm;
    AdjointGradient ag;

};

int main(int argc, char *argv[]) {

    #if MT
        int th = 4;
        omp_set_num_threads(th);
        Eigen::setNbThreads(th);
    #endif

    double H = 1.0;
    double W = 1.0;
    int VW = 32;
    int VH = 32;
    double Q = 20.0/0.001;
    double Cmet = 65.0;
    double Cpla = 0.2;
    double M = 0.4;
    int p = 1;
    int p_max = 5;

    #if TIME
    int max_iter = 1;
    #else
    int max_iter = 250;
    #endif

    std::vector<double> v(VW*VH, 0.0), v_old(v);

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
    double f, g, f_old = 0;
    std::vector<double> df(n), dg(n);
    HeatEq heateq(H, W, VW, VH, Q, Cmet, Cpla, p, M, BcBottom, BcRight, BcTop, BcLeft);

    int m = 1;
    double movlim = 0.5;
    double V_MIN = 0.0, V_MAX = 1.0;
    std::vector<double> v_min(n, V_MIN), v_max(n, V_MAX);
    MMASolver *mma = new MMASolver(n,m);
    double ch = 1.0;

	int itr = 0;
	while (itr < max_iter) {
		itr++;
		heateq(v, f, df, g, dg);

		// Set outer move limits
		for (int i=0;i<n;i++) {
			v_max[i] = min(V_MAX, v[i] + movlim);
			v_min[i] = max(V_MIN, v[i] - movlim);
        }

        #if TIME
		    auto t_start = std::chrono::system_clock::now();
        #endif

        mma->Update(v.data(), df.data(), &g, dg.data(), v_min.data(), v_max.data());

        #if TIME
            auto t_end = std::chrono::system_clock::now();
            std::chrono::duration<double> diff = t_end-t_start;
            std::cout << "MMA up took: " << diff.count() << " s" << std::endl;
        #endif

		ch = inf_norm_diff(v, v_old);
		v_old = v;

        // Update after x iterations. Better for large grids.
//		if (itr % 50 == 0 && p < p_max){
//		    heateq.update_p(++p);
//		}

        // Update if objective is not changing and constraint is below threshold
        if (abs((f - f_old)) / n < 1e-6 && g < 1e-2 && p < p_max){
            heateq.update_p(++p);
        }

        f_old = f;
		printf("it.: %d, obj.: %f, g.: %f, ch.: %f \n",itr, f, g, ch);
	}


	std::vector<double> t(n);
	heateq.solve_heat_eq(v, t);

	std::ofstream file;
	file.open("../results/result.out");

	// Write v to file
	for (int i = 0; i < n -1; i++) {
	    file << v[i] << ",";
	}
	file << v[n - 1] << "\n";

	// Write t to file
    for (int i = 0; i < n -1; i++) {
        file << t[i] << ",";
    }
    file << t[n - 1] << "\n";

	file.close();

	// Deallocate
	delete mma;
	return 0;

}


