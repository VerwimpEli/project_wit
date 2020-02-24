clear;clc;

%%Bestand met enkele van de test 
 VB = 24;VH = 24; % Aantal volumes in de hoogte en breedte. Incluisief de kleinere op de randen
 Varray = ones(VB,VH); %Volledig metaal materiaal
 Varray(1/3*VB:2/3*VB-1,:) = zeros(VB/3,VH);
 Q = 0;
 BC0 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Onder geisoleerde rand
 BC1 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]]; % Rechter 
 BC2 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Boven geisoleerde rand
 BC3 = [['D',1,1,0];['D',2,VH-1,0];['D',VH,VH,0]];% Linker 
[Sol,K] = FVM(VB,VH,Varray,Q,65,0.2,BC0,BC1,BC2,BC3);
 SOL = reshape(Sol,[VB,VH]);
 figure(1); surf(SOL); 
 title("Testcase 1 Stroken verschillend materiaal");

%%%%TESTCASE 2 %Domein 1mx1m
VH = 5; VB = 24; % Aantal volumes in de hoogte en breedte. Incluisief de kleinere op de randen
Varray = ones(VB,VH); %Volledig metaal materiaal
Q = 200;
BC0 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Onder geisoleerde rand
BC1 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]]; % Rechter 
BC2 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Boven geisoleerde rand
BC3 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]];% Linker 


[Sol,K] = FVM(VB,VH,Varray,Q,65,0.2,BC0,BC1,BC2,BC3);
SOL = reshape(Sol,[VB,VH]);
figure(2); surf(SOL);

Sol_Theo = T_theo(20,20,Q,1,24,65);

figure(3);
plot(SOL(:,1)-Sol_Theo');

%%TESTCASE 2B nieuwe


%%%%TEST CASE 3 --> Convergentie van TEST CASE 2
%Mesh = [3,4,5,6,7,8,9,10,11,12,16,24,48,96,192]; %VB
Mesh = 5:5:250;
Dx = 1./(Mesh-1);
Error = zeros(size(Mesh));

for i = 1:size(Mesh,2)
    VB = Mesh(i);VH = 10;
    Varray = ones(Mesh(i),10); %Volledig metaal materiaal
    BC0 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Onder geisoleerde rand
    BC1 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]]; % Rechter 
    BC2 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Boven geisoleerde rand
    BC3 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]];% Linker 

    [Sol,K] = FVM(Mesh(i),VH,Varray,Q,65,0.2,BC0,BC1,BC2,BC3);
    SOL = reshape(Sol,[Mesh(i),VH]);
    Sol_Theo = T_theo(0,0,Q,1,Mesh(i),65);
    %Error(i) = max(abs(SOL(:,1)-Sol_Theo')); MaxNorm
    Error(i) = norm(abs(SOL(:,1)-Sol_Theo'))/Mesh(i); %grid 2-Norm
end
%semilogy(Dx,Error); %verkeerde soort plot voor het juiste verband te zien
figure(4); loglog(Dx,Error); hold on; grid on; 
title("2De orde convergentie voor paraboolisch testprobleem");
xlabel("\Deltax"); ylabel("2-Norm van de error / n");

function Sol_Theo = T_theo(T1,T2,Q,L,VH,K)
    x = linspace(0,L,VH);
    Sol_Theo = -Q.*x.*x/2/K + (T2-T1+Q*L^2/2/K)/L.*x + T1;
end