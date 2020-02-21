clear;
%%ProtoType van de Solver voor het oplossen van de warmte vergelijkingen.
%%De gebruikte is Finite elements

%Stap 1 Inlezen van de Mesh en eigenschappen van het materiaal.
load('Node-Element_QUAD.mat'); %Laden van de reeds aangemaakte Node en element arrays 
load('GaussN2.mat'); % 
NoNode = size(NodeArray,1); %elke node heeft 1 DOF
NoElement = size(ElementArray,1);
%De node array heeft volgende structuur NodeNummer - Xcoord - Ycoord
%De element array heeft volgende structuur ElementNummer - Node1 - Node2 -
%Node3 - Materiaalverhouding(1 metaal). De ordening van de elementen is
%in wijzerzin

%Constanten
Cmet = 65; %Metaal
Cpla = 0.2; %Plastiek


%%%%%%%%%%%%%%%%%%%%%%%%
%Stap2 Opstellen van de stijfheidsmatrix
K = zeros(NoNode,NoNode);
%Dit gebeurt element per element. De evaluatie van de integraal gebeurt
%door middel van gauss kwadratuur voor driehoek (of eventueel vierkant voor
%latere elementen

%Voor vierkante elementen %
for i = 1:NoElement
    C = Cmet*ElementArray(i,6)+Cpla*(1-ElementArray(i,6)); %De warmte coefficient in een element
    C = C/0.001; %De dikte in rekening brengen
    
    KE = zeros(4);
    for j = 1:size(Gauss,1)
        %Elementstijfheidsmatrix berekenen ADHV Gauss midpoint rule
        R = calcR(Gauss(j,1),Gauss(j,2));
        CO = [[NodeArray(ElementArray(i,2),2),NodeArray(ElementArray(i,2),3)];[NodeArray(ElementArray(i,3),2),NodeArray(ElementArray(i,3),3)];[NodeArray(ElementArray(i,4),2),NodeArray(ElementArray(i,4),3)];[NodeArray(ElementArray(i,5),2),NodeArray(ElementArray(i,5),3)]];
        J = R*CO; %De matrix J (deze zou normaal diagonaal zijn voor de gestructureerd mesh)
    
        AM = J\R;
        KE = KE + Gauss(j,3)*C*transpose(AM)*AM*abs(det(J));
    end
    
    
    %De coefficienten bijtellen is de globale matrix
    Index = [ElementArray(i,2),ElementArray(i,3),ElementArray(i,4),ElementArray(i,5)];
    K(Index,Index) = K(Index,Index) + KE;
end

%spy(K);

%%%%%%%%%%%%%%%%%%%%%%%%
%Stap3 Invullen van "Right hand Side" van het stelsel
RHS = zeros(NoNode,1);
% volumetrisch/oppervlakte warmteproductie 
q = 2/(0.01*0.01);
for i = 1:NoElement
    %Niet efficient, tegelijk bereken met stijfheids matrix
    R = calcR(0,0);
    CO = [[NodeArray(ElementArray(i,2),2),NodeArray(ElementArray(i,2),3)];[NodeArray(ElementArray(i,3),2),NodeArray(ElementArray(i,3),3)];[NodeArray(ElementArray(i,4),2),NodeArray(ElementArray(i,4),3)];[NodeArray(ElementArray(i,5),2),NodeArray(ElementArray(i,5),3)]];
    J = R*CO; %De matrix J (deze zou normaal diagonaal zijn voor de gestructureerd mesh)
    
    RHS(ElementArray(i,2:5)) = RHS(ElementArray(i,2:5)) + 4*q*abs(det(J))*[1/4;1/4;1/4;1/4]; %Midpoint rule %lijkt ok
end

%%%%%%%%%%%%%%%%%%%%%%%%
%Stap5 implementeren van de Dirichlet randvoorwaarden
T = 293; 
%index_c = [1,5,12,16]; u_c = [T-20;T-20;T;T];
index_c = [12,16]; u_c = [T;T];
%Stap5B Implementeren van de Neumann voorwaarden
%Homogene --> moet eigelijks niets veranderen

%Stap5C Schrappen van de vergelijkingen
index_f = setdiff( 1:NoNode, index_c);
K_ff = K(index_f,index_f);
K_fc = K(index_f,index_c);

%%%%%%%%%%%%%%%%%%%%%%%%
%Stap6 Oplossen vh het systeem
u_f = K_ff\(RHS(index_f)-K_fc*u_c);

%%%%%%%%%%%%%%%%%%%%%%%%
%Stap7 Visualisatie van resultaten.
Sol = zeros(NoNode,1);
Sol(index_c) = u_c; Sol(index_f) = u_f;

%%%%%%%%%%%%%%%%%%%%%%%%
%Plotten van Quads
figure;hold on;
xlabel('X'); ylabel('Y'); zlabel('Temperatuur');
for i = 1:NoElement
    %4 lijnen per element plotten 
    for j = 2:4
        N1 = ElementArray(i,j); N2 = ElementArray(i,j+1);
        plot3([NodeArray(N1,2),NodeArray(N2,2)],[NodeArray(N1,3),NodeArray(N2,3)],[Sol(N1),Sol(N2)],'k');
    end
    
    N1 = ElementArray(i,2); N2 = ElementArray(i,5);
    plot3([NodeArray(N1,2),NodeArray(N2,2)],[NodeArray(N1,3),NodeArray(N2,3)],[Sol(N1),Sol(N2)],'k');
end

%%Functies
function R = calcR(xksi,eta) %vanuit oefenzitten NMM
    R = 1/4*[[-(1-eta);-(1-xksi)],[-(1+eta);(1-xksi)],[(1+eta);(1+xksi)],[(1-eta);-(1+xksi)]]; %Mogelijks niet in de juiste volgorde
end