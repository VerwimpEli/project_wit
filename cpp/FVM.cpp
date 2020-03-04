//////////////////////////////////////////////
///Finite volume methode
//////////////////////////////////////////////

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

// INPUTS voor functie 
K :: Referentie naar spaarse matrix van V
SOL :: referentie naar solution vector
*/

//////OUTPUTS NA

//Defines en includes
#define NEUMANN 1
#define DIRICHLET 0
#include <vector> //Voor de vectoren te kunnen gebruiken
#include <math.h>       

//FUNCTOR
template<typename S> //type scalar dat gebruikt wordt voor berekeneningen
class FVM
{
    public:
    FVM(float H, float W, int VW, int VH, float M, float Q, float Cmet, float Cpla, int p,)
    :H_(H);W_(W);VW_(VW);VH_(VH);M_(M);Q_(Q); Cmet_(Cmet); Cpla_(Cpla); p_(p);
    diag_(VW*VH);diagU1_(VW*VH);diagU2_(VW*VH);RHS_(VW*VH,Q*W_/(VW-1)*H_/(VH-1));dx_(W_/(VW-1));dy_(H_/(VH-1));
    {}//Constructor

    void operator()(std::vector<double> & v, Eigen::SparseMatrix<?> & K, std::vector<double> & SOL)
    {
        //Hulpvariable
        S C;
        S M;
        std::vector<S> Material(v); //Zou een kopie moeten zijn

        //Berekenen van de materiaal coefficienten //MatArray = (1 - v) .^ p * Cpla + v .^ p * Cmet;
        std::transform(Material.begin(),Material.end(),Material.begin(),[float Cmet, float Cpla](double ve) -> std::pow((1-ve),p_)*Cpla_ + std::pow(ve,p_)*Cmet_ ))

        //Berekenen van de RHS //Verschalen randen links en rechts want dit zijn kleinere volumes
        for(int k = 0;k<VW*VH; k = k + VW) 
        {
            RHS(k) = RHS(k)/2;
            RHS(k+VW-1) = RHS(k+VW-1)/2;
        }

        for(int k = 0; k<VW;k = k + 1) //Verschalen randen onder en boven 
        {
            RHS(k) = RHS(k)/2;
            RHS(k+VW*(VH-1) = RHS(k+VW*(VH-1))/2;
        }

        //Contributies voor de horizontale interface tussen cellen
        for(int i =0;i<VW-1 ;i = i+1) //Over de breedte
        {
            

        }

    }



private: 
//Diagonalen
std::vector<S> diag_;
std::vector<S> diagU1_;
std::vector<S> diagU2_;
std::vector<S> RHS_;
//Overige parameters
float H_;
float W_;
int VW_;
int VH_;
float M_;
float Q_;
float Cmet_;
float Cpla_;
int p_;
S dx_;
S dy_;
};

