//////////////////////////////////////////////
///Finite volume methode     /////////////////
//////////////////////////////////////////////
#ifndef CPP_FVM
#define CPP_FVM

//Defines en includes
#include <assert.h> //Controle ect
#include <iostream> //Input en output
#include <vector> //Voor de vectoren te kunnen gebruiken
#include "Eigen/Sparse" //Solver
#include "Eigen/SparseQR"
#include "Eigen/Dense"
#include "Eigen/OrderingMethods"
#include <math.h>  // std::pow
#include <cmath>  // std::floor
#include <algorithm> //Transform
#include "util.h"
#include "BoundaryCondition.h"
#ifdef UmfLUSolver
#include <Eigen/UmfPackSupport>
#endif

/**
 * FVM functor: calculates the temperature with the finite volume method on a square grid. The material constants
 * are calculated using a local harmoinc average of v, to prevent checkerboards. A p value is also used to steer
 * the solutions towards 0's and 1's. The constructor takes the problem specific parameter. There are two () operators,
 * one will only calculate a solution, the other will also calculate an L vector, which is required by the adjoint
 * gradient method. This happens here because factorizing of K happens here as well.
 *
 * The functors require:
 * - v: the current solution
 * - Sol: A place to store the solution
 * - K: A solution to store the K matrix
 * optional:
 * - L: a place to store the solution of KL = -1 for the adjoint.
 *
 */

template<typename S> //type scalar dat gebruikt wordt voor berekeneningen
class FVM
{
    private: 

    S const H_; //Hoogte vh vierkant domein
    S const W_; //Breedte vh vierkant domein
    int const VW_; //Aantal cellen in de breedte
    int const VH_; //aantal cellen in de hoogte
    S const Q_; //oppervlakte warmte productie
    S const Cmet_; //materiaal coefficient van metaal
    S const Cpla_; // materiaal coefficient van plasiek
    int p_;
    S const dx_;
    S const dy_;

    // Diagonalen
    std::vector<S> diag_; //De hoofd diagonaal
    std::vector<S> diagU1_; // de eerste, nevendiagonaal
    std::vector<S> diagU2_; // de verder gelegen diagonaal op afstand VW_ gelegen
    std::vector<S> RHS_; // Rechter zijde vh stelsel

    // De boundary objecten
    BoundaryCondition const BC0_;
    BoundaryCondition const BC1_;
    BoundaryCondition const BC2_;
    BoundaryCondition const BC3_;

    // functie om voor een functorcall alle gebruikte vector (voor de constructie van K) te reseten
    void reset(){
        std::fill(diag_.begin(), diag_.end(), 0.0);
        std::fill(diagU1_.begin(), diagU1_.end(), 0.0);
        std::fill(diagU2_.begin(), diagU2_.end(), 0.0);
    }
    
    public:
    FVM(S H, S W, int VW, int VH, S Q, S Cmet, S Cpla, int p, BoundaryCondition BC0, BoundaryCondition BC1, BoundaryCondition BC2, BoundaryCondition BC3)
    :H_(H), W_(W), VW_(VW), VH_(VH), Q_(Q), Cmet_(Cmet), Cpla_(Cpla), p_(p), dx_(W_/(VW-1)), dy_(H_/(VH-1)),
    diag_(VW*VH), diagU1_(VW*VH), diagU2_(VW*VH), RHS_(VW*VH,Q*W_/(VW-1)*H_/(VH-1)),
    BC0_(BC0), BC1_(BC1), BC2_(BC2), BC3_(BC3)
    {

    }//Constructor


    void operator()(std::vector<double> const & v,  std::vector<double> & SOL, Eigen::SparseMatrix<double, 0> & K){
        Eigen::VectorXd L;
        (*this)(v, SOL, K, L, false);
    }

    void operator()(std::vector<double> const & v,  std::vector<double> & SOL, Eigen::SparseMatrix<double, 0> & K,
            Eigen::VectorXd & L, bool l_solve = true)
    {
        reset();

        //Hulpvariable
        S C;
        S M;
        std::vector<S> Material(v); //Zou een kopie moeten zijn

        //Berekenen van de materiaal coefficienten //MatArray = (1 - v) .^ p * Cpla + v .^ p * Cmet;
        std::transform(Material.begin(),Material.end(),Material.begin(),
                [Cmet_c = Cmet_, Cpla_c = Cpla_,p_c = p_](S ve) ->
                S {return std::pow((1-ve),p_c)*Cpla_c + std::pow(ve,p_c)*Cmet_c;});


        std::fill(RHS_.begin(), RHS_.end(), Q_*dx_*dy_); //Initialisatie van RHS
        //Berekenen van de RHS //Verschalen randen links en rechts want dit zijn kleinere volumes
        for(int k = 0; k<VW_*VH_; k = k + VW_)
        {
            RHS_[k] = RHS_[k]/2;
            RHS_[k+VW_-1] = RHS_[k+VW_-1]/2;
        }

        for(int k = 0; k<VW_;k = k + 1) //Verschalen randen onder en boven
        {
            RHS_[k] = RHS_[k]/2;
            RHS_[k+VW_*(VH_-1)] = RHS_[k+VW_*(VH_-1)]/2;
        }


        //Contributies voor de horizontale interface tussen cellen
        for(int i =0;i<VW_-1 ;i = i+1) //Over de breedte
        {
            //Cell gelegen op de onderrand
            M = 2/(1/Material[i] + 1/Material[i+1]);
            C = M*dy_/2/dx_;
            diag_[i] = diag_[i] + C; // k  = i //let op, dit is de matlab notatie
            diag_[i+1] = diag_[i+1] + C; // k = i+1 //let op, dit is de matlab notatie
            diagU1_[i] = diagU1_[i] - C; // k = i //let op, dit is de matlab notatie

            //Interne cell
            for(int j = 1; j<VH_-1;j = j +1)
            {
                M = 2/(1/Material[i+j*VW_] + 1/Material[i+1+j*VW_]);
                C = M*dy_/dx_;
                diag_[i + j*VW_] = diag_[i + j*VW_] + C; // k  = i+VW*(j-1) //let op, dit is de matlab notatie
                diag_[i + j*VW_+1] = diag_[i + j*VW_+1] + C; // k = i+VW*(j-1)+1 //let op, dit is de matlab notatie
                diagU1_[i + j*VW_] = diagU1_[i + j*VW_] - C; // k = i+VW*(j-1) //let op, dit is de matlab notatie
            }
            
            //Cellen gelegen op de boven rand
            M = 2/(1/Material[i+(VH_-1)*VW_] + 1/Material[i+1+(VH_-1)*VW_]);
            C = M*dy_/2/dx_;
            diag_[i+(VH_-1)*VW_] = diag_[i+(VH_-1)*VW_] + C; // k  = i+VW*(VH-1) //let op, dit is de matlab notatie
            diag_[i+(VH_-1)*VW_+1] = diag_[i+(VH_-1)*VW_+1] + C; // k = i+VW*(VH-1)+1 //let op, dit is de matlab notatie
            diagU1_[i+(VH_-1)*VW_] = diagU1_[i+(VH_-1)*VW_] - C; // k = i+VW*(VH-1) //let op, dit is de matlab notatie
        }


        //Contributies voor de verticale interface tussen de cellen
        for(int j = 0; j<VH_-1;j = j +1)
        {
            //Cellen op de Linkerzijde
            M = 2/(1/Material[j*VW_] + 1/Material[(j+1)*VW_]);
            C = M*dx_/2/dy_;

            diag_[j*VW_] = diag_[j*VW_] + C; // k = 1+VW*(j-1); //let op, dit is de matlab notatie
            diag_[(j+1)*VW_] = diag_[(j+1)*VW_] + C; 
            diagU2_[j*VW_] = diagU2_[j*VW_] - C; // k = 1+VW*(j-1); //let op, dit is de matlab notatie

            //Interne Cellen
            for(int i = 1; i < VW_-1; i = i +1)
            {
                M = 2/(1/Material[i+j*VW_] + 1/Material[i+(j+1)*VW_]);
                C = M*dx_/dy_;

                diag_[i+j*VW_] = diag_[i+j*VW_] + C; // k = i+VW*(j-1);; //let op, dit is de matlab notatie
                diag_[i+(j+1)*VW_] = diag_[i+(j+1)*VW_] + C; 
                diagU2_[i+j*VW_] = diagU2_[i+j*VW_] - C; // k = i+VW*(j-1);; //let op, dit is de matlab notatie
            }

            //cellen op de rechterzijde
            M = 2/(1/Material[VW_-1+j*VW_] + 1/Material[VW_-1+(j+1)*VW_]);
            C = M*dx_/2/dy_;

            diag_[VW_-1+j*VW_] = diag_[VW_-1+j*VW_] + C; // VW+VW*(j-1); //let op, dit is de matlab notatie
            diag_[VW_-1+(j+1)*VW_] = diag_[VW_-1+(j+1)*VW_] + C; 
            diagU2_[VW_-1+j*VW_] = diagU2_[VW_-1+j*VW_] - C; // VW+VW*(j-1); //let op, dit is de matlab notatie
            
        }
        


        //AFHANDELEN VD BOUNDARY CONDITIONS
        S PW = std::pow(10,7)*(*max_element(diag_.begin(), diag_.end()));//Penaltywaarde voor de DIRCHLET randvoorwaarden te implementeren met de penalty methode
        int beginIndex;
        int endIndex;
        //BC0 - Onders
        //De kleine cell in het begin vd rand
        if(BC0_.GetStart().type() == DIRICHLET)
        {
            RHS_[0] = RHS_[0] + BC0_.GetStart().value()*PW; //k = 1 //matlabnotatie
            diag_[0] = diag_[0] + PW;
        }
        else
        {
            RHS_[0] = RHS_[0] + BC0_.GetStart().value()*dx_/2;
        }
        //De Gewone cellen op de onderrand
        for (BoundarySegment seg: BC0_.GetSegments())
        {
            
            beginIndex = std::max(1,ratioToIndex(seg.start(),VW_));
            endIndex = std::min(VW_-1,ratioToIndex(seg.stop(),VW_));
            if(seg.type() == DIRICHLET)
            {
                
                for( int i = beginIndex; i < endIndex;i = i + 1)
                {
                    RHS_[i] = RHS_[i] + seg.value()*PW; //k = i //matlabnotatie
                    diag_[i] = diag_[i] + PW; 
                }
                
            }
            else
            {
                for( int i = beginIndex; i < endIndex;i = i + 1)
                {
                    RHS_[i] = RHS_[i] + seg.value()*dx_;
                }
            }
        }
        //De kleine cell op het einde vd rand
        if(BC0_.GetStop().type() == DIRICHLET)
        {
            RHS_[VW_] = RHS_[VW_] + BC0_.GetStop().value()*PW; //k = VW //matlabnotatie
            diag_[VW_] = diag_[VW_] + PW;
        }
        else
        {
            RHS_[VW_] = RHS_[VW_] + BC0_.GetStop().value()*dx_/2;
        }

        //Test
        //Print(diag_); //Ziet er oke uit voor homogene Neumann (dit zegt eigenlijk niet veel)
        //Print(RHS_);  //Ziet er oke uit voor homogene Neumann

        //BC1 - RECHTS
        //De kleine cell in het begin vd rand
        if(BC1_.GetStart().type() == DIRICHLET)
        {
            RHS_[VW_-1] = RHS_[VW_-1] + BC1_.GetStart().value()*PW; //k = 1 //matlabnotatie
            diag_[VW_-1] = diag_[VW_-1] + PW;
        }
        else
        {
            RHS_[VW_-1] = RHS_[VW_-1] + BC1_.GetStart().value()*dy_/2;
        }
         //De Gewone cellen op de rechter
        for (BoundarySegment seg: BC1_.GetSegments())
        {
            
            beginIndex = std::max(1,ratioToIndex(seg.start(),VH_));
            endIndex = std::min(VH_-1,ratioToIndex(seg.stop(),VH_));
            if(seg.type() == DIRICHLET)
            {
                
                for( int i = beginIndex; i < endIndex;i = i + 1)
                {
                    RHS_[VW_-1 + i*VW_] = RHS_[VW_-1 + i*VW_] + seg.value()*PW; 
                    diag_[VW_-1 + i*VW_] = diag_[VW_-1 + i*VW_] + PW; 
                }
                
            }
            else
            {
                for( int i = beginIndex; i < endIndex;i = i + 1)
                {
                    RHS_[VW_-1 + i*VW_] = RHS_[VW_-1 + i*VW_] + seg.value()*dy_;
                }
            }
        }
        //Kleine cell op het einde vd rand
        if(BC1_.GetStop().type() == DIRICHLET)
        {
            RHS_[VH_*VW_-1] = RHS_[VH_*VW_-1] + BC1_.GetStop().value()*PW; 
            diag_[VH_*VW_-1] = diag_[VH_*VW_-1] + PW;
        }
        else
        {
            RHS_[VH_*VW_-1] = RHS_[VH_*VW_-1] + BC1_.GetStop().value()*dy_/2;
        }

        //BC2 - Boven
        //De kleine cell in het begin vd rand
        if(BC2_.GetStart().type() == DIRICHLET)
        {
            RHS_[VW_*(VH_-1)] = RHS_[VW_*(VH_-1)] + BC2_.GetStart().value()*PW; //k = 1 //matlabnotatie
            diag_[VW_*(VH_-1)] = diag_[VW_*(VH_-1)] + PW;
        }
        else
        {
            RHS_[VW_*(VH_-1)] = RHS_[VW_*(VH_-1)] + BC2_.GetStart().value()*dx_/2;
        }
         //De Gewone cellen op de onderrand
        for (BoundarySegment seg: BC2_.GetSegments())
        {
            
            beginIndex = std::max(1,ratioToIndex(seg.start(),VW_));
            endIndex = std::min(VW_-1,ratioToIndex(seg.stop(),VW_));
            if(seg.type() == DIRICHLET)
            {
                
                for( int i = beginIndex; i < endIndex;i = i + 1)
                {
                    RHS_[i + VW_*(VH_-1)] = RHS_[i + VW_*(VH_-1)] + seg.value()*PW; 
                    diag_[i + VW_*(VH_-1)] = diag_[i + VW_*(VH_-1)] + PW; 
                }
                
            }
            else
            {
                for( int i = beginIndex; i < endIndex;i = i + 1)
                {
                    RHS_[i + VW_*(VH_-1)] = RHS_[i + VW_*(VH_-1)] + seg.value()*dx_;
                }
            }
        }
        //Kleine cell op het einde vd rand
        if(BC2_.GetStop().type() == DIRICHLET)
        {
            RHS_[VH_*VW_-1] = RHS_[VH_*VW_-1] + BC2_.GetStop().value()*PW; 
            diag_[VH_*VW_-1] = diag_[VH_*VW_-1] + PW;
        }
        else
        {
            RHS_[VH_*VW_-1] = RHS_[VH_*VW_-1] + BC2_.GetStop().value()*dx_/2;
        }


        //BC3 - LINKS
        //De kleine cell in het begin vd rand
        if(BC3_.GetStart().type() == DIRICHLET)
        {
            RHS_[0] = RHS_[0] + BC3_.GetStart().value()*PW; //k = 1 //matlabnotatie
            diag_[0] = diag_[0] + PW;
        }
        else
        {
            RHS_[0] = RHS_[0] + BC3_.GetStart().value()*dy_/2;
        }
         //De Gewone cellen op de onderrand
        for (BoundarySegment seg: BC3_.GetSegments())
        {

            beginIndex = std::max(1,ratioToIndex(seg.start(),VH_));
            endIndex = std::min(VH_-1,ratioToIndex(seg.stop(),VH_));
            if(seg.type() == DIRICHLET)
            {

                for( int i = beginIndex; i < endIndex;i = i + 1)
                {
                    RHS_[i*VW_] = RHS_[i*VW_] + seg.value()*PW;
                    diag_[i*VW_] = diag_[i*VW_] + PW;
                }

            }
            else
            {
                for( int i = beginIndex; i < endIndex;i = i + 1)
                {
                    RHS_[i*VW_] = RHS_[i*VW_] + seg.value()*dy_;
                }
            }
        }
        //Kleine cell op het einde vd rand
        if(BC3_.GetStop().type() == DIRICHLET)
        {
            RHS_[(VH_-1)*VW_] = RHS_[(VH_-1)*VW_] + BC3_.GetStop().value()*PW;
            diag_[(VH_-1)*VW_] = diag_[(VH_-1)*VW_] + PW;
        }
        else
        {
            RHS_[(VH_-1)*VW_] = RHS_[(VH_-1)*VW_] + BC3_.GetStop().value()*dy_/2;
        }
        
        // Constructie van K uit de diagonalen
        K.reserve(Eigen::VectorXi::Constant(VH_*VW_, 5));
        for (int i = 0; i < VW_ * VH_; i++){
            K.insert(i, i) = diag_[i];
            if (i < VW_ * VH_ - 1) {
                K.insert(i, i+1) = diagU1_[i];
                K.insert(i + 1, i) = diagU1_[i];
            }
            if (i < VW_ * VH_ - VW_) {
                K.insert(i, i + VW_) = diagU2_[i];
                K.insert(i + VW_, i) = diagU2_[i];
            }
        }

        K.makeCompressed();
        Eigen::VectorXd RHS_E(VW_ * VH_);
        Eigen::VectorXd x;

        for (int i  = 0; i < VW_ * VH_; i++){
            RHS_E(i) = RHS_[i];
        }

        // VW <= 256: SimplicialLDLT is fast enough
        // VW >= 256: Use UmfPackLU. (SuiteSparse etc required). Has multi threading support.

        #ifdef UmfLUSolver
        Eigen::UmfPackLU<Eigen::SparseMatrix<double>> solver;
        #else
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        #endif

        solver.compute(K);
        x = solver.solve(RHS_E);

        if (l_solve){
            Eigen::VectorXd rhs(VW_ * VH_);
            for (int i = 0; i < VW_ * VH_; i++) {
                rhs(i) = -1;
            }
            scale(rhs, VW_, rhs.rows());
            L = solver.solve(rhs);
        }

        for (int i = 0; i < VW_ * VH_; i++){
            SOL[i] = x(i);
        }
    }

    void update_p(int p){
        p_ = p;
    }
};

#endif
