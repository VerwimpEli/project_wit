clear;  close all;

%Script dat enkele testen uitvoerd om de correctheid van de functie
%Harmonic_FVM te testen.
%Script for verification of the code.
%In order to run this make sure matlab has acces to the other files aswell
%by "adding to path" the tests folder 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                TEST CASE 1 : 3 STRIPS OF MATERIAL       %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3 verschillende stroken materiaal ter controlen dat de de juiste materiaal
%parameters geselecteerd worden. 
%Simulation of plate consiting of 3 stripes of metal. There should be 3 distict regions with the
%gradient over the plastic region. This verifies material selection.
 VB = 24;VH = 24; % Aantal volumes in de hoogte en breedte. Incluisief de kleinere op de randen
 H = 1; B = 1; %Hoogte en breedte van het domein
 Varray = ones(VB,VH); %Volledig metaal materiaal
 Varray(1/3*VB:2/3*VB-1,:) = zeros(VB/3,VH);
 Q = 0;
 BC0 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Onder geisoleerde rand
 BC1 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]]; % Rechter 
 BC2 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Boven geisoleerde rand
 BC3 = [['D',1,1,0];['D',2,VH-1,0];['D',VH,VH,0]];% Linker 
[Sol,K] = Harmonic_FVM(B,H,VB,VH,Varray,Q,65,0.2,BC0,BC1,BC2,BC3);
 SOL = reshape(Sol,[VB,VH]);
 figure(1); surf(SOL); 
 title("Testcase 1 Harmonic_FVM Stroken verschillend materiaal");
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                TEST CASE 2 : Parabool                  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulatie van een uniforme plaat met aan 2 zijde isolatie en 2 zijde een
%opgelegde temperatuur. Door inwendige warmte generatie wordt er een
%parabool profiel aangemaakt door FVM. De berekende waarden liggen exact op
%dit profiel
%Simulation of a uniform plate with on 2 opposide side isolation (neumann =
%0) and 2 given temperatures. As result of the uniform heat production
%within the plate, the solution is a parabola. 
VH = 5; VB = 24; % Aantal volumes in de hoogte en breedte. Incluisief de kleinere op de randen
H = 1; B = 2; %Hoogte en breedte van het domein 
Varray = ones(VB,VH); %Volledig metaal materiaal
Q = 200;
BC0 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Onder geisoleerde rand
BC1 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]]; % Rechter 
BC2 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Boven geisoleerde rand
BC3 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]];% Linker 


[Sol,K] = Harmonic_FVM(B,H,VB,VH,Varray,Q,65,0.2,BC0,BC1,BC2,BC3);
SOL = reshape(Sol,[VB,VH]);
figure(2); surf(SOL); title("Test Case 2 : Parabool gegenereerd door Harmonic_FVM");

Sol_Theo = T_theo(20,20,Q,B,24,65); %Vergelijken met theortisch 
figure(3);
plot((SOL(:,1)-Sol_Theo')./Sol_Theo'); title("Test Case 2 : verschil tussen  Harmonic_FVM en Theoretisch profiel");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                TEST CASE 3 : Convergentie van Test Case 2 %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %De 2de testcase wordt herhaald met verschillende meshes. Gezien de
 %eigenschappen van de methode zou de parabool telkens exact geinterpoleerd
 %moeten worden %dit is dus gewoon nul berekenen (Met numerieke errors)!
 %1ste orde Convergentie gedrag is het gevolg van dynamische penalty waarde
 %(Lineair met mesh/ proportioneel met grootste waarde op diagonaal
 %Test case 2 for different meshes. Result should still be zero for all. 
Mesh = 5:5:1500;
Dx = 1./(Mesh-1);
Error = zeros(size(Mesh));

for i = 1:size(Mesh,2)
    VB = Mesh(i);VH = 10;
    Varray = ones(Mesh(i),10); %Volledig metaal materiaal
    BC0 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Onder geisoleerde rand
    BC1 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]]; % Rechter 
    BC2 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Boven geisoleerde rand
    BC3 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]];% Linker 

    [Sol,K] = Harmonic_FVM(B,H,Mesh(i),VH,Varray,Q,65,0.2,BC0,BC1,BC2,BC3);
    SOL = reshape(Sol,[Mesh(i),VH]);
    Sol_Theo = T_theo(20,20,Q,B,Mesh(i),65);
    Error(i) = max(abs((SOL(:,1)-Sol_Theo')./Sol_Theo')); %MaxNorm
    %Error(i) = norm(abs(SOL(:,1)-Sol_Theo'))/Mesh(i); %grid 2-Norm
end
figure(4); loglog(Dx,Error); hold on; grid on; 
title("Test Case 2 voor verschillende mesh afstanden");
xlabel("\Deltax"); ylabel("Max-Norm van de REL error ");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TEST CASE 4 : NEUMANN : BC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Testcase for a non-zero neumann Boundary condition. The resulting
%solutions is linear. And the solutions found are therefor exact. relative error should 
%be zero (Not numerically of course)
VH = 5; VB = 24; % Aantal volumes in de hoogte en breedte. Incluisief de kleinere op de randen
H = 1; B = 2; %Hoogte en breedte van het domein 
Varray = ones(VB,VH); %Volledig metaal materiaal
Q = 0;
BC0 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Onder geisoleerde rand
BC1 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]]; % Rechter 
BC2 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Boven geisoleerde rand
BC3 = [['N',1,1,20];['N',2,VH-1,20];['N',VH,VH,20]];% Linker 


[Sol,K] = Harmonic_FVM(B,H,VB,VH,Varray,Q,65,0.2,BC0,BC1,BC2,BC3);
SOL = reshape(Sol,[VB,VH]);
figure(5); surf(SOL); title("TestCase 4 : Voorbeeld van een oplossing");

Sol_Theo = T_theo_neumann(-20/65,20,Q,B,24,65);

figure(6);
plot((SOL(:,1)-Sol_Theo')./Sol_Theo');  title("TestCase 4 : Relative Error met theoretische oplossing");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TEST CASE 5 : NEUMANN : BC - Convergentie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Testcase 4 for different meshes. The relative error should always be zero(Numerically).
%1ste orde convergence behavior is the consquence of the dynamische penalty value
%(Lineair with mesh/ proportional with largest value of the diagonal
Mesh = 5:5:1000;
Dx = 1./(Mesh-1);
Error = zeros(size(Mesh));

for i = 1:size(Mesh,2)
    VB = Mesh(i);VH = 10;
    Varray = ones(Mesh(i),10); %Volledig metaal materiaal
    BC0 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Onder geisoleerde rand
    BC1 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]]; % Rechter 
    BC2 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Boven geisoleerde rand
    BC3 = [['N',1,1,20];['N',2,VH-1,20];['N',VH,VH,20]];% Linker 

    [Sol,K] = Harmonic_FVM(B,H,Mesh(i),VH,Varray,Q,65,0.2,BC0,BC1,BC2,BC3);
    SOL = reshape(Sol,[Mesh(i),VH]);
    Sol_Theo = T_theo_neumann(-20/65,20,Q,B,Mesh(i),65);
    Error(i) = max(abs((SOL(:,1)-Sol_Theo')./Sol_Theo')); %MaxNorm
    
end
%semilogy(Dx,Error); %verkeerde soort plot voor het juiste verband te zien
figure(7); loglog(Dx,Error); hold on; grid on; 
title("TestCase 5 : convergentie voor testprobleem met neumann randvoorwaarde");
xlabel("\Deltax"); ylabel("Max-Norm van de REL error ");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TEST CASE 6 : Niet uniforme Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convergentie for a Uniform material with non uniform heat production
%The exact solution is a 3de orde polynomial. Due to the non unformness of
%the heat production the solution is no longer exact. 2nd order
%convergece should be seen.
VH = 5; H = 1; B = 1;
Mesh = 5:10:2000;
Dx = 1./(Mesh-1);
Error = zeros(size(Mesh));

for i = 1:size(Mesh,2)
    VB = Mesh(i);
    Varray = ones(Mesh(i),5); %Volledig metaal materiaal
    Q = CreateForcingMatrix_Q(VB,VH,65,1); %aanmaken van de forcing matrix
    BC0 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Onder geisoleerde rand
    BC1 = [['D',1,1,293];['D',2,VH-1,293];['D',VH,VH,293]]; % Rechter 
    BC2 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Boven geisoleerde rand
    BC3 = [['D',1,1,273];['D',2,VH-1,273];['D',VH,VH,273]];% Linker 

    [Sol,K] = Harmonic_FVM_NonUniformQ(B,H,Mesh(i),VH,Varray,Q,65,0.2,BC0,BC1,BC2,BC3);
    SOL = reshape(Sol,[Mesh(i),VH]);
    Sol_Theo = T_theo_poly(Mesh(i),65,1);
    Error(i) = max(abs((SOL(:,1)-Sol_Theo')./Sol_Theo')); %MaxNorm
    %Error(i) = norm(abs(SOL(:,1)-Sol_Theo'))/Mesh(i); %grid 2-Norm
end
%semilogy(Dx,Error); %verkeerde soort plot voor het juiste verband te zien
figure(8); loglog(Dx,Error); hold on; grid on; 
title({'2nd order convergence of the relative error of the FVM temperature',' for a 4th order polynomial testproblem'});
xlabel("\Deltax"); ylabel("Max-Norm of the REL error");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%HULP FUNCTIES %%%%%%%%%%%%%%%%%

function Sol_Theo = T_theo(T1,T2,Q,L,VH,K) %Test Parabool
    x = linspace(0,L,VH);
    Sol_Theo = -Q.*x.*x/2/K + (T2-T1+Q*L^2/2/K)/L.*x + T1;
end

function Sol_Theo = T_theo_neumann(DT1,T2,Q,L,VH,K) %Test Neumann
    x = linspace(0,L,VH);
    Sol_Theo = -Q.*x.*x/2/K + DT1.*x -DT1*L + T2 + Q/2/K*L^2;
end

% %aanmaken van de matrix Q
% function Q = CreateForcingMatrix_Q(VB,VH,K,L)
%     dx = L/(VB-1);
%     Q = zeros(VB,VH);
%     for i = 1:VB
%             Q(i,:) = -90*K*(dx/2 + (i-1)*dx);
%     end
% end

%aanmaken van de matrix Q
function Q = CreateForcingMatrix_Q(VB,VH,K,L)
    dx = L/(VB-1);
    Q = zeros(VB,VH);
    for i = 2:VB-1
            Q(i,:) = -240*K*((i-1)*dx).^2;
    end
    
     %Begin Rand
     Q(1,:) = -240*K*(0.25*dx).^2;
     %Eind rand
     Q(VB,:) = -240*K*(L-0.25*dx).^2;
end

%Theoretische functie
function Sol_Theo = T_theo_poly(V,K,L)
    x = linspace(0,L,V);
    Sol_Theo = 20.*x.*x.*x.*x  + 273;
end