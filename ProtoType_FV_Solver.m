clear;

%Constanten %Mesh wordt volledig vierkant verondersteld
H = 0.005; B = 0.005; %Hoogte en breedte van het domein
VH = 25; VB = 25; % Aantal volumes in de hoogte en breedte
dx = B/VB; dy = H/VH; %Cell grotes
T = 293; %Temperatuur Dirichlet
q = 2;
Cmet = 65; %Metaal
Cpla = 0.2; %Plastiek

%MateriaalArray
MatArray = Cmet*ones(VB,VH);

%Aanmaken van matrix en RHS
K = zeros(VB*VH);
RHS = ones(VB*VH,1);

%temperatuursproductie in de cell
RHS= q*dx*dy*RHS;

%Stensil voor de inwendige punten %Geen efficient implementatie (sommige
%resultaten kunnen voor meerder )
%Zowel temperatuur en materiaal state zijn cell centered
%Verticale Faces
for i = 1:VB-1 %in breedte
   for j = 1:VH %in hooghte
       C1 = MatArray(i,j)*dx/dy; C2 = MatArray(i+1,j)*dx/dy;
       k = i+VB*(j-1);
       %Vergelijking k
       K(k,k) = K(k,k) - C1; K(k,k+1) = K(k,k+1) + C2;
       %vergelijking k+1
       K(k+1,k) = K(k+1,k) + C1; K(k+1,k+1) = K(k+1,k+1) - C2;
   end
end

%Horizontale Faces
for i = 1:VB %in breedte
   for j = 1:VH-1 %in hooghte
       C1 = MatArray(i,j)*dx/dy; C2 = MatArray(i,j+1)*dx/dy;
       k = i+VB*(j-1);
       %Vergelijking k
       K(k,k) = K(k,k) - C1; K(k,k+VB) = K(k,k+VB) + C2;
       %vergelijking k+VB
       K(k+VB,k) = K(k+VB,k) + C1; K(k+VB,k+VB) = K(k+VB,k+VB) - C2;
   end
end

%%%%%%%%%%%%%%%ùboundaryConditions
%Neumann -> niets voor homogene boundaryconditions
%Dirichlet -> 'Backwardsdifference'

%Linkse Rand
T1 = 0;
for i =1:VH
    C1 = 2*MatArray(1,i)*dx/dy; C2 = 2*T1*MatArray(1,i)*dx/dy;
    k = 1+(i-1)*VB;
    K(k,k) = K(k,k) - C1; RHS(k) = RHS(k) - C2;
end

%Rechtse rand 
T2 = 0;
for i =1:VH
    C1 = 2*MatArray(VB,i)*dx/dy; C2 = 2*T2*MatArray(VB,i)*dx/dy;
    k = i*VB;
    K(k,k) = K(k,k) - C1; RHS(k) = RHS(k) - C2;
end

%Oplossen Systeem
Sol = K\RHS;

%Visualisatie
Sol = reshape(Sol,[VB,VH]);
surf(Sol);


