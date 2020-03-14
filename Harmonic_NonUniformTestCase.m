%Constanten %Mesh wordt volledig vierkant verondersteld %Speciale
%Vierkanten aan de randen zodat de berekende temperaturen het volledige
%domein insluiten
%Mesh

Mesh = 5:5:500;
Dx = 1./(Mesh-1);
Error = zeros(size(Mesh));

for i = 1:size(Mesh,2)
    VB = Mesh(i);VH = 5;
    Varray = ones(Mesh(i),5); %Volledig metaal materiaal
    Q = CreateForcingMatrix_Q(VB,VH,65,1); %aanmaken van de forcing matrix
    BC0 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Onder geisoleerde rand
    BC1 = [['D',1,1,288];['D',2,VH-1,288];['D',VH,VH,288]]; % Rechter 
    BC2 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Boven geisoleerde rand
    BC3 = [['D',1,1,273];['D',2,VH-1,273];['D',VH,VH,273]];% Linker 

    [Sol,K] = Harmonic_FVM_NonUniformQ(Mesh(i),VH,Varray,Q,65,0.2,BC0,BC1,BC2,BC3);
    SOL = reshape(Sol,[Mesh(i),VH]);
    Sol_Theo = T_theo(Mesh(i),65,1);
    Error(i) = max(abs(SOL(:,1)-Sol_Theo')); %MaxNorm
    %Error(i) = norm(abs(SOL(:,1)-Sol_Theo'))/Mesh(i); %grid 2-Norm
end
%semilogy(Dx,Error); %verkeerde soort plot voor het juiste verband te zien
figure(4); loglog(Dx,Error); hold on; grid on; 
title("1Ste orde convergentie voor paraboolisch testprobleem");
xlabel("\Deltax"); ylabel("Max-Norm van de error / n");

%aanmaken van de matrix Q
function Q = CreateForcingMatrix_Q(VB,VH,K,L)
    dx = L/(VB-1);
    Q = zeros(VB,VH);
    for i = 1:VB
            Q(i,:) = -90*K*(dx/2 + (i-1)*dx);
    end
end
%Theoretische functie
function Sol_Theo = T_theo(V,K,L)
    x = linspace(0,L,V);
    Sol_Theo = 15.*x.*x.*x  + 273;
end