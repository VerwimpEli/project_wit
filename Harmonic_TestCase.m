clear;

%Constanten %Mesh wordt volledig vierkant verondersteld %Speciale
%Vierkanten aan de randen zodat de berekende temperaturen het volledige
%domein insluiten
VH = 20; VB = 20; % Aantal volumes in de hoogte en breedte. Incluisief de kleinere op de randen
%rng(500); %Reproduceerbaarheid
Varray = rand(VB*VH,1);
v = reshape(Varray, [VB, VH]);
%%%Boundary condition
BC0 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Onder geisoleerde rand
BC1 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]]; % Rechter 
BC2 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Boven geisoleerde rand
BC3 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]];% Linker 

%%Extra parameters
q = 200;
Cmet = 65; Cpla = 0.2;
p = 2;

%Refentie oplossing
[Sol,K] = Harmonic_FVM(VB, VH, v, q, Cmet, Cpla, BC0, BC1, BC2, BC3, p);


%Visualisatie %Heel ruw komt niet direct overeen met echte systeem
SOL = reshape(Sol,[VB,VH]);
figure(1);
surf(SOL); 
xlabel("X"); ylabel("Y"); zlabel("Temperatuur")

%Benaderen Jacobiaan
 FDG = FD_G(VB,VH,Varray,q,Cmet, Cpla, BC0, BC1, BC2, BC3, p)';

%Adjoint
L = (K')\-scale(ones(VB*VH,1));
AG = Harmonic_Adjoint_Gradient(VB,VH,v,L,Sol, p, Cmet, Cpla);

ERR1 = log10(abs(reshape(AG-FDG,[VB,VH])./reshape(AG+FDG,[VB,VH])/2));


figure(2);
subplot(1,3,1); 
surf(reshape(AG,[VB,VH])); 
title("AG");
subplot(1,3,2)
surf(ERR1); 
title("REL ERR1 TOV AVG");
xlabel("X"); ylabel("Y"); zlabel("ERROR1")
subplot(1,3,3)
surf(reshape(FDG,[VB,VH])); 
title("FDG");

%%%Bereken van de gradient DMV Finite difference
function J = FD_G(VB,VH,Varray,q,Cmet, Cpla, BC0, BC1, BC2, BC3, p)
    Delta = 10^-6;
    J = zeros(size(Varray));
    %Referentie Oplossing
    [Sol, ~, ~, ~] =  Harmonic_heateq(Varray, 1, VB, VH, q, Cmet, Cpla, BC0, BC1, BC2, BC3, p);
    
    %Loop
    for i = 1:size(Varray,1)
            Varray2 = Varray; Varray2(i) = Varray2(i)*(1+Delta);
            [FD_Sol, ~, ~, ~] =  Harmonic_heateq(Varray2, 1, VB, VH, q, Cmet, Cpla, BC0, BC1, BC2, BC3, p);
            J(i)= (FD_Sol - Sol)/(Delta*Varray(i)); %of is dit de gradient?
    end
end
