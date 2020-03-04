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

//FUNCTOR
template<typename S> //type scalar dat gebruikt wordt
class FVM
{
    public:
    FVM(float H, float W, int VW, int VH, float M, float Q, float Cmet, float Cpla, int p,)
    :H_(H);W_(W);VW_(VW);VH_(VH);M_(M);Q_(Q); Cmet_(Cmet); Cpla_(Cpla); p_(p);
    diag_(VW*VH);diagU1_(VW*VH);diagL1_(VW*VH);diagU2_(VW*VH);diagL2_(VW*VH);
    {}//Constructor

    void operator()(std::vector<double> & v, Eigen::SparseMatrix<?> & K, std::vector<double>)
    {
        S dx = W_/(VW-1);
        S dy = H_/(VH-1);
    }



private: //Hulpvariable
//Diagonalen
std::vector<S> diag_;
std::vector<S> diagU1_;
std::vector<S> diagL1_;
std::vector<S> diagU2_;
std::vector<S> diagL2_;
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
};

