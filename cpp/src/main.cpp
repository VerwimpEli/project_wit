#include <chrono>
#include <iostream>
#include "MMASolver.h"
#include "util.h"
#include "BoundaryCondition.h"
#include "adjoint.cpp"
#include "FVM.cpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <fstream>
#include <numeric>
#include <iterator>
#include <vector>
#include <algorithm>

/**
 * This functor calculates all parameters required for the optimizer.
 * Given v, it will calculate:
 *  - f: the value of the objective function as sum(T). See FVM.cpp.
 *  - df: the gradient of the objective. See adjoint.cpp.
 *  - g: the value of the constraint: sum(v) - M*n^2.
 *  - dg: the gradient of the constraint
 */
class HeatEq
{
public:
    HeatEq(double H, double W, int VW, int VH, double Q, double Cmet, double Cpla, int p, double M,
           BoundaryCondition BC0, BoundaryCondition BC1, BoundaryCondition BC2, BoundaryCondition BC3)
            : H_(H), W_(W), Cmet_(Cmet), Cpla_(Cpla), Q_(Q), M_(M),  VW_(VW), VH_(VH), p_(p),
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

        #ifdef TIME
            auto t_end = std::chrono::system_clock::now();
            std::chrono::duration<double> diff = t_end-t_start;
            std::cout << "FVM took: " << diff.count() << " s" << std::endl;
        #endif

        t_start = std::chrono::system_clock::now();
        ag(v, L, sol, df);

        #ifdef TIME
            t_end = std::chrono::system_clock::now();
            diff = t_end-t_start;
            std::cout << "ADJ took: " << diff.count() << " s" << std::endl;
        #endif

        scale(sol, VW_, sol.size());
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

    auto t_start_full = std::chrono::system_clock::now();

    int max_iter, VW, VH;

    if(argc == 1){
        max_iter = 500;
        VW = 32;
        VH = 32;
    } else if (argc == 2){
        max_iter = std::atoi(argv[1]);
        VW = 32;
        VH = 32;
    } else if (argc == 3){
        max_iter = std::atoi(argv[1]);
        VW = std::atoi(argv[2]);
        VH = VW;
    } else {
        std::cout << "Too many arguments" << std::endl;
        exit(1);
    }

    double H = 1.0;
    double W = 1.0;
    double Q = 2/0.001;
    double Cmet = 65.0;
    double Cpla = 0.2;
    double M = 0.4;
    int p = 1;
    int p_max = 5;

    #ifdef TIME
    max_iter = 1;
    #endif

//
    std::vector<double> v(VW*VH, 0.0), v_old(v);

//    Read and interpolate old v

//    std::fstream is("../results/result_106_32.out", std::ios::in);
//    std::istream_iterator<double> start(is), end;
//    std::vector<double> v_read(start, end);
//    std::cout << v_read.size() << std::endl;
//    std::vector<double> v = interpolate(v_read, VW/2), v_old(v);
//    std::vector<double> v(v_read), v_old(v);

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

	int itr = 0;
	while (itr < max_iter) {
		itr++;
		heateq(v, f, df, g, dg);

		// Set outer move limits
		for (int i=0;i<n;i++) {
			v_max[i] = min(V_MAX, v[i] + movlim);
			v_min[i] = max(V_MIN, v[i] - movlim);
        }

        #ifdef TIME
        auto t_start = std::chrono::system_clock::now();
        #endif

        mma->Update(v.data(), df.data(), &g, dg.data(), v_min.data(), v_max.data());

        #ifdef TIME
        auto t_end = std::chrono::system_clock::now();
        std::chrono::duration<double> diff = t_end-t_start;
        std::cout << "MMA took: " << diff.count() << " s" << std::endl;
        #endif

		v_old = v;

        // Update if objective is not changing and constraint is below threshold
        if (abs((f - f_old)) / n < 1e-4 && g < 1e-2 && p < p_max){
            heateq.update_p(++p);
        }

        f_old = f;
#ifndef TIME
	    heateq(v, f, dg, g, dg);

    printf("it.: %d, obj.: %f, g.: %f, ch.: %f \n",itr, f, g, inf_norm_diff(v, v_old));
#endif
	}

    auto t_end_full = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = t_end_full-t_start_full;
    std::cout << "Program took: " << diff.count() << " s" << std::endl;


#ifndef TIME

    // Write v and t to file
	std::vector<double> t(n);
	heateq.solve_heat_eq(v, t);

	std::ofstream file;
    std::time_t tm = std::time(0);
    std::tm* now = std::localtime(&tm);
    char f_name[40];

    std::strftime(f_name, 40, "../results/result_%j_%H%M%S.out", now);

	file.open(f_name);

	// Write v to file
	for (int i = 0; i < n -1; i++) {
	    file << v[i]  << ",";
	}
	file << v[n - 1] << "\n";

	// Write t to file
    for (int i = 0; i < n -1; i++) {
        file << t[i] << ",";
    }
    file << t[n - 1] << "\n";

	file.close();

#endif
	// Deallocate
	delete mma;
	return 0;

}


