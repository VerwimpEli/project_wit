#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <math.h>
#include <vector>
#include <valarray>
#include <functional> 
#include "BoundaryCondition.cpp"

class AdjointGradient {

    private:
        float W_;
        float H_;
        int VW_;
        int VH_;
        float M_;
        float Q_;
        float Cmet_;
        float Cpla_;
        int p_;

    public: 
        AdjointGradient(float W, float H, int VW, int VH, float M, float Q, float Cmet,
        float Cpla, int p)
            : W_(W)
            , H_(H)
            , VW_(VW)
            , VH_(VH)
            , M_(M)
            , Q_(Q)
            , Cmet_(Cmet)
            , Cpla_(Cpla)
            , p_(p)
        {};

#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCDFAInspection"
        void operator() (std::vector<double> v, std::vector<double> L,
        std::vector<double> SOL, std::vector<double> AG){
            float dx = W_/(VW_-1);
            float dy = H_/(VH_-1);
            // for_each(v.begin(), v.end(), [](double i) 
            // { 
            //     1.0-i; 
            // });
            // std::transform(v.begin(), v.end(), v.begin(),
            //    std::bind(std::minus<T>(), std::placeholders::_1, 1));
            // std::vector<double> MatArray = std::pow(v, 2)*;

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

            //Berekenen van de materiaal coefficienten
            std::transform(MatArray.begin(),MatArray.end(),MatArray.begin(), [Cmet_c = Cmet_, Cpla_c = Cpla_,p_c = p_](double ve)
                    -> double {return std::pow((1-ve),p_c)*Cpla_c + std::pow(ve,p_c)*Cmet_c;});

            std::transform(MatDerArray.begin(),MatDerArray.end(),MatDerArray.begin(), [Cmet_c = Cmet_, Cpla_c = Cpla_,p_c = p_](double ve)
                    -> double {return -p_c*std::pow((1-ve),p_c-1)*Cpla_c + p_c*std::pow(ve,p_c-1)*Cmet_c;});

            for(int i = 0; i<VW_-1; ++i){
                int k = i;
                M1 = MatArray[i*VW_]; M2 = MatArray[(i+1)*VW_];
                dM1 = 2/(pow(M1,2)*pow(1/M1+1/M2,2)); dM2 = 2/(pow(M2,2)*pow(1/M1+1/M2,2));

                R1 = MatDerArray[i*VW_]*dM1*dy/dx/2; R2 = MatDerArray[(i+1)*VW_]*dM2*dy/dx/2;

                dg[k] = dg[k] + R1*SOL[k] - R1*SOL[k+1];
                dg[k+1] = dg[k+1] - R2*SOL[k] + R2*SOL[k+1];
                subdg[k] = subdg[k] - R1*SOL[k] + R1*SOL[k+1];
                supdg[k+1] = supdg[k+1] + R2*SOL[k] - R2*SOL[k+1];

                for(int j=1; j<VH_-1; ++j){
                    k = i+VW_*(j-1);

                    M1 = MatArray[i*VW_+j]; M2 = MatArray[(i+1)*VW_+j];
                    dM1 = 2/(pow(M1,2)*pow(1/M1+1/M2,2)); dM2 = 2/(pow(M2,2)*pow(1/M1+1/M2,2));

                    R1 = MatDerArray[i*VW_+j]*dM1*dy/dx/2; R2 = MatDerArray[(i+1)*VW_+j]*dM2*dy/dx/2;

                    dg[k] = dg[k] + R1*SOL[k] - R1*SOL[k+1];
                    dg[k+1] = dg[k+1] - R2*SOL[k] + R2*SOL[k+1];
                    subdg[k] = subdg[k] - R1*SOL[k] + R1*SOL[k+1];
                    supdg[k+1] = supdg[k+1] + R2*SOL[k] - R2*SOL[k+1];
                }

                k = i+VW_*(VH_-1);

                M1 = MatArray[i*VW_+VH_]; M2 = MatArray[(i+1)*VW_+VH_];
                dM1 = 2/(pow(M1,2)*pow(1/M1+1/M2,2)); dM2 = 2/(pow(M2,2)*pow(1/M1+1/M2,2));

                R1 = MatDerArray[i*VW_+VH_]*dM1*dy/dx/2; R2 = MatDerArray[(i+1)*VW_+VH_]*dM2*dy/dx/2;

                dg[k] = dg[k] + R1*SOL[k] - R1*SOL[k+1];
                dg[k+1] = dg[k+1] - R2*SOL[k] + R2*SOL[k+1];
                subdg[k] = subdg[k] - R1*SOL[k] + R1*SOL[k+1];
                supdg[k+1] = supdg[k+1] + R2*SOL[k] - R2*SOL[k+1];
            }

            for(int i = 0; i<VH_-1; ++i){
                int k = 1+VW_*(i-1);

                M1 = MatArray[i]; M2 = MatArray[i+1];
                dM1 = 2/(pow(M1,2)*pow(1/M1+1/M2,2)); dM2 = 2/(pow(M2,2)*pow(1/M1+1/M2,2));

                R1 = MatDerArray[i]*dM1*dy/dx/2; R2 = MatDerArray[i+1]*dM2*dy/dx/2;

                dg[k] = dg[k] + R1*SOL[k] - R1*SOL[k+VW_];
                dg[k+VW_] = dg[k+VW_] - R2*SOL[k] + R2*SOL[k+VW_];
                subvwdg[k] = subvwdg[k] - R1*SOL[k] + R1*SOL[k+VW_];
                supvwdg[k+VW_] = supvwdg[k+VW_] + R2*SOL[k] - R2*SOL[k+VW_];

                for(int j=1; j<VW_-1; ++j){
                    k = j+VW_*(i-1);

                    M1 = MatArray[j*VW_+i]; M2 = MatArray[j*VW_+i+1];
                    dM1 = 2/(pow(M1,2)*pow(1/M1+1/M2,2)); dM2 = 2/(pow(M2,2)*pow(1/M1+1/M2,2));

                    R1 = MatDerArray[j*VW_+i]*dM1*dy/dx/2; R2 = MatDerArray[j*VW_+i+1]*dM2*dy/dx/2;

                    dg[k] = dg[k] + R1*SOL[k] - R1*SOL[k+VW_];
                    dg[k+VW_] = dg[k+VW_] - R2*SOL[k] + R2*SOL[k+VW_];
                    subvwdg[k] = subvwdg[k] - R1*SOL[k] + R1*SOL[k+VW_];
                    supvwdg[k+VW_] = supvwdg[k+VW_] + R2*SOL[k] - R2*SOL[k+VW_];
                }

                k = VW_+VW_*(i-1);

                M1 = MatArray[VW_*(VW_-1)+i]; M2 = MatArray[i+1+VW_*(VW_-1)];
                dM1 = 2/(pow(M1,2)*pow(1/M1+1/M2,2)); dM2 = 2/(pow(M2,2)*pow(1/M1+1/M2,2));

                R1 = MatDerArray[VW_*(VW_-1)+i]*dM1*dy/dx/2; R2 = MatDerArray[i+1+VW_*(VW_-1)]*dM2*dy/dx/2;

                dg[k] = dg[k] + R1*SOL[k] - R1*SOL[k+VW_];
                dg[k+1] = dg[k+1] - R2*SOL[k] + R2*SOL[k+VW_];
                subvwdg[k] = subvwdg[k] - R1*SOL[k] + R1*SOL[k+VW_];
                supvwdg[k+VW_] = supvwdg[k+VW_] + R2*SOL[k] - R2*SOL[k+VW_];
            }

        }
#pragma clang diagnostic pop


};