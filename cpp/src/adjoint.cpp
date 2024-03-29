#ifndef CPP_ADJOINT
#define CPP_ADJOINT

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <math.h>
#include <vector>
#include <valarray>
#include <functional>
#include "Eigen/Sparse"

/**
 * This functor will calculate the adjoint gradient for the heat problem.
 * The constructor takes the problem specific parameters. The () operator needs:
 *  - v: the current solutions
 *  - L: solution of KL = -1, calculated in FVM
 *  - Sol: the current temperature solution calculated by FVM
 *  - AG: a vector to store the result
 */
class AdjointGradient {

    private:
        float H_;
        float W_;
        int VW_;
        int VH_;
        float M_;
        float Q_;
        float Cmet_;
        float Cpla_;
        int p_;

    public: 
        AdjointGradient(float H, float W, int VW, int VH, float M, float Q, float Cmet,
        float Cpla, int p)
            : H_(H)
            , W_(W)
            , VW_(VW)
            , VH_(VH)
            , M_(M)
            , Q_(Q)
            , Cmet_(Cmet)
            , Cpla_(Cpla)
            , p_(p)
        {};

        void operator() (std::vector<double> const & v, Eigen::VectorXd & L, std::vector<double> & SOL, std::vector<double> & AG){

            float dx = W_/(VW_-1);
            float dy = H_/(VH_-1);

            std::vector<double> MatArray(v);
            std::vector<double> MatDerArray(v);
            std::vector<double> dg(VH_*VW_, 0);
            std::vector<double> subdg(VH_*VW_, 0);
            std::vector<double> supdg(VH_*VW_, 0);
            std::vector<double> subvwdg(VH_*VW_, 0);
            std::vector<double> supvwdg(VH_*VW_, 0);
            double M1;
            double M2;
            double dM1;
            double dM2;
            double R1;
            double R2;

            std::transform(MatArray.begin(),MatArray.end(),MatArray.begin(), [Cmet_c = Cmet_, Cpla_c = Cpla_,p_c = p_](double ve)
                    -> double {return std::pow((1-ve),p_c)*Cpla_c + std::pow(ve,p_c)*Cmet_c;});

            std::transform(MatDerArray.begin(),MatDerArray.end(),MatDerArray.begin(), [Cmet_c = Cmet_, Cpla_c = Cpla_,p_c = p_](double ve)
                    -> double {return -p_c*std::pow((1-ve),p_c-1)*Cpla_c + p_c*std::pow(ve,p_c-1)*Cmet_c;});


            for(int i = 0; i<VW_-1; ++i){

                int k = i;
                M1 = MatArray[i];
                M2 = MatArray[i+1];

                dM1 = 2/(pow(M1,2)*pow(1/M1+1/M2,2));
                dM2 = 2/(pow(M2,2)*pow(1/M1+1/M2,2));

                R1 = MatDerArray[i] * dM1 * dy/ dx /2;
                R2 = MatDerArray[i+1] * dM2* dy / dx /2;

                dg[k] = dg[k] + R1*SOL[k] - R1*SOL[k+1];
                dg[k+1] = dg[k+1] - R2*SOL[k] + R2*SOL[k+1];
                subdg[k] = subdg[k] - R1*SOL[k] + R1*SOL[k+1];
                supdg[k+1] = supdg[k+1] + R2*SOL[k] - R2*SOL[k+1];


                for(int j=1; j<VH_-1; ++j){

                    k = i+VW_*j;

                    M1 = MatArray[i + VW_*j]; M2 = MatArray[(i+1) + VW_*j];
                    dM1 = 2/(pow(M1,2)*pow(1/M1+1/M2,2));
                    dM2 = 2/(pow(M2,2)*pow(1/M1+1/M2,2));

                    R1 = MatDerArray[i + VW_*j]*dM1*dy/dx;
                    R2 = MatDerArray[(i+1) + VW_*j]*dM2*dy/dx;

                    dg[k] = dg[k] + R1*SOL[k] - R1*SOL[k+1];
                    dg[k+1] = dg[k+1] - R2*SOL[k] + R2*SOL[k+1];
                    subdg[k] = subdg[k] - R1*SOL[k] + R1*SOL[k+1];
                    supdg[k+1] = supdg[k+1] + R2*SOL[k] - R2*SOL[k+1];
                }

                k = i+VW_*(VH_-1);

                M1 = MatArray[i + VW_*(VH_-1)]; M2 = MatArray[(i+1) + VW_*(VH_-1)];
                dM1 = 2/(pow(M1,2)*pow(1/M1+1/M2,2));
                dM2 = 2/(pow(M2,2)*pow(1/M1+1/M2,2));

                R1 = MatDerArray[i+ VW_*(VH_-1)]*dM1*dy/dx/2;
                R2 = MatDerArray[(i+1)+VW_*(VH_-1)]*dM2*dy/dx/2;

                dg[k] = dg[k] + R1*SOL[k] - R1*SOL[k+1];
                dg[k+1] = dg[k+1] - R2*SOL[k] + R2*SOL[k+1];
                subdg[k] = subdg[k] - R1*SOL[k] + R1*SOL[k+1];
                supdg[k+1] = supdg[k+1] + R2*SOL[k] - R2*SOL[k+1];
            }

            for(int i = 0; i<VH_-1; ++i){

                int k = VW_*i;

                M1 = MatArray[i*VH_]; M2 = MatArray[(i+1) * VH_];
                dM1 = 2/(pow(M1,2)*pow(1/M1+1/M2,2));
                dM2 = 2/(pow(M2,2)*pow(1/M1+1/M2,2));

                R1 = MatDerArray[i * VH_]*dM1*dy/dx/2;
                R2 = MatDerArray[(i+1) * VH_]*dM2*dy/dx/2;


                dg[k] = dg[k] + R1*SOL[k] - R1*SOL[k + VW_];
                dg[k + VW_] = dg[k + VW_] - R2*SOL[k] + R2*SOL[k + VW_];
                subvwdg[k] = subvwdg[k] - R1*SOL[k] + R1*SOL[k + VW_];
                supvwdg[k + VW_] = supvwdg[k + VW_] + R2*SOL[k] - R2*SOL[k + VW_];


                for(int j=1; j<VW_-1; ++j){
                    k = j+VW_*i;

                    M1 = MatArray[j + (VW_*i)];
                    M2 = MatArray[j + VW_*(i+1)];
                    dM1 = 2/(pow(M1,2)*pow(1/M1+1/M2,2));
                    dM2 = 2/(pow(M2,2)*pow(1/M1+1/M2,2));

                    R1 = MatDerArray[j + VW_*i]*dM1*dy/dx;
                    R2 = MatDerArray[j + VW_*(i+1)]*dM2*dy/dx;

                    dg[k] = dg[k] + R1*SOL[k] - R1*SOL[k + VW_];
                    dg[k + VW_] = dg[k + VW_] - R2*SOL[k] + R2*SOL[k + VW_];
                    subvwdg[k] = subvwdg[k] - R1*SOL[k] + R1*SOL[k + VW_];
                    supvwdg[k + VW_] = supvwdg[k + VW_] + R2*SOL[k] - R2*SOL[k + VW_];
                }

                k = VW_+VW_*i-1;

                M1 = MatArray[(VW_-1)+ VW_* i];
                M2 = MatArray[(VW_-1)+ VW_*(i+1)];
                dM1 = 2/(pow(M1,2)*pow(1/M1+1/M2,2)); dM2 = 2/(pow(M2,2)*pow(1/M1+1/M2,2));

                R1 = MatDerArray[(VW_-1)+VW_*i]*dM1*dy/dx/2;
                R2 = MatDerArray[(VW_-1)+VW_*(i+1)]*dM2*dy/dx/2;

                dg[k] = dg[k] + R1*SOL[k] - R1*SOL[k + VW_];
                dg[k + VW_] = dg[k + VW_] - R2*SOL[k] + R2*SOL[k + VW_];
                subvwdg[k] = subvwdg[k] - R1*SOL[k] + R1*SOL[k + VW_];
                supvwdg[k + VW_] = supvwdg[k + VW_] + R2*SOL[k] - R2*SOL[k + VW_];
            }

            Eigen::SparseMatrix<double> G(VW_ * VH_, VW_ * VH_);
            G.reserve(Eigen::VectorXi::Constant(VH_*VW_, 5));

            for (int i = 0; i < VH_ * VW_; i++){
                G.insert(i, i) = dg[i];
                if (i < VW_ * VH_ - 1) {
                    G.insert(i, i+1) = supdg[i+1];
                    G.insert(i + 1, i) = subdg[i];
                }
                if (i < VW_ * VH_ - VW_) {
                    G.insert(i, i + VW_) = supvwdg[i+VW_];
                    G.insert(i + VW_, i) = subvwdg[i];
                }
            }

            Eigen::VectorXd adjgrad(VW_ * VH_);
            adjgrad = L.transpose() * G;

            for (int i = 0; i < VW_ * VH_; i++){
                AG[i] = adjgrad(i);
            }
        }

        void update_p(int p){
            p_ = p;
        }
};

#endif
