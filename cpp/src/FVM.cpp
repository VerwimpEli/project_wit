//////////////////////////////////////////////
///Finite volume methode     /////////////////
//////////////////////////////////////////////
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

#ifndef TIME
#define TIME 0;
#endif

#ifndef DIRICHLET
#define  DIRICHLET 0
#endif

#ifndef NEUMANN
#define  NEUMANN 1
#endif


//VERVOLG FVM

////BESCHRIJVING
/*  Functor die de finite volume methode geassocieerd met een vierkant domein oplost. De mesh is gebaseerd op gebruikte mesh 
in de paper. De mesh voor de state en design variabelen zijn het zelfde. Het harmonische gemiddelde is gebruikt voor het 
berekenen van de materiaal constante op de cell interfaces, omdat het checkerboard pateroon hierdoor vermeden wordt. De mesh 
begint bij element 1 in de linker onderhoek, element 2 bevind zich rechts, ect. Op het einde van de onderste rij, wordt naar
de rechterzijde van de op één na onderste rij gesprongen.

*/

///////INPUTS voor initialsatie
/* 
H, W :: de Hoogte (H) en Widt (W) van het domein
VH, VW :: het aantal volumes gebruikt voor het discritiseren van de hoogte (VH) en de breedte (VW)
v :: std-vector met daarin dit materiaal ratios van 2 materialen in een element, in dezelfde ordening als de mesh
q :: de uniforme warmte productie van het oppervlak gegeven in w/m²
Cmet :: coefficient van het metaal/materiaal1 (de 1 in de v)
Cpla :: coefficient van het plastiek/materiaal2 (de 0 in de v)  
BC0,BC1,BC2,BC3 --> bijvoorbeeld uit matlab
BC0 = [['N',1,1,0];['N',2,VW-1,0];['N',VW,VW,0]]; %Onder geisoleerde rand
BC1 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]]; % Rechter 
BC2 = [['N',1,1,0];['N',2,VW-1,0];['N',VW,VW,0]]; %Boven geisoleerde rand
BC3 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]];% Linker 
p :: de simp parameter

// INPUTS voor functie /functor
v :: 
K :: Referentie naar spaarse matrix van V
SOL :: referentie naar solution vector
*/

//////OUTPUTS NA

//Functie voor het bepalen van een positie op de rand (tussen 0 en 1) naar de index vh element.
//De functie houdt rekening met de kleine cell op de rand
int ratioToIndex(double ratio, int Cells) //Pass by value
{
    double dx = 1.0/(Cells-1);
    if(ratio <= dx/2){ return 0;}
    else if( ratio >= 1-dx/2){return Cells -1;}
    else{ return std::floor((ratio-dx/2)/dx) + 1;}
}

//FUNCTOR
template<typename S> //type scalar dat gebruikt wordt voor berekeneningen
class FVM
{
    private: 
    //Diagonalen
    std::vector<S> diag_; //De hoofd diagonaal
    std::vector<S> diagU1_; // de eerste, nevendiagonaal
    std::vector<S> diagU2_; // de verder gelegen diagonaal op afstand VW_ gelegen
    std::vector<S> RHS_; // Rechter zijde vh stelsel

    //Overige parameters
    int p_;

    S const H_; //Hoogte vh vierkant domein
    S const W_; //Breedte vh vierkant domein
    int const VW_; //Aantal cellen in de breedte
    int const VH_; //aantal cellen in de hoogte
    S const Q_; //oppervlakte warmte productie
    S const Cmet_; //materiaal coefficient van metaal
    S const Cpla_; // materiaal coefficient van plasiek
    S const dx_;
    S const dy_;
    BoundaryCondition const BC0_; //De boundary objecten
    BoundaryCondition const BC1_;
    BoundaryCondition const BC2_;
    BoundaryCondition const BC3_;

    //functie om voor een functorcall alle gebruikte vector (voor de constructie van K) te reseten
    void reset(){
        std::fill(diag_.begin(), diag_.end(), 0.0);
        std::fill(diagU1_.begin(), diagU1_.end(), 0.0);
        std::fill(diagU2_.begin(), diagU2_.end(), 0.0);
    }
    
    public:
    FVM(S H, S W, int VW, int VH, S Q, S Cmet, S Cpla, int p, BoundaryCondition BC0, BoundaryCondition BC1, BoundaryCondition BC2, BoundaryCondition BC3)
    :H_(H), W_(W), VW_(VW), VH_(VH), Q_(Q), Cmet_(Cmet), Cpla_(Cpla), p_(p),
    diag_(VW*VH), diagU1_(VW*VH), diagU2_(VW*VH), RHS_(VW*VH,Q*W_/(VW-1)*H_/(VH-1)), dx_(W_/(VW-1)), dy_(H_/(VH-1)),
    BC0_(BC0), BC1_(BC1), BC2_(BC2), BC3_(BC3)
    {
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
        //std::cout<<"RHS"<<std::endl;
        //Print(RHS_); //Ziet er goed uit
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
        S PW = std::pow(10,8); //Penaltywaarde voor de DIRCHLET randvoorwaarden te implementeren met de penalty methode
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

//        std::cout << "\n------------ FVM ---------------" << std::endl;
//        Print(diag_);
//        Print(diagU1_);
//        Print(diagU2_);
//        std::cout << "\n---------- FVM end -------------" << std::endl;

        //Constructie van K uit de diagonalen
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

        #if UmfLUSolver
            Eigen::UmfPackLU<Eigen::SparseMatrix<double>> solver;
        #else
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        #endif

        #if TIME
            auto t_start = std::chrono::system_clock::now();
        #endif
        solver.compute(K);

        #if TIME
            auto t_end = std::chrono::system_clock::now();
            std::chrono::duration<double> diff = t_end-t_start;
            std::cout << "Compute took: " << diff.count() << " s" << std::endl;
        #endif

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
