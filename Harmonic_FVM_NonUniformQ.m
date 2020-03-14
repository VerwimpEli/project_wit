%De functie FVM, deze functie is de simulatie van de differentiaal
%vergelijking zelf. De mesh is in overeenkomst met de paper over heat
%topolisation/optimasation. elk volume heeft zijn eigen
%materiaal(coefficent)
%Inputs
%VW --> Int,Aantal volumes in de breedte/width richting
%VH --> Int,aantal volumes in de hoogste
%v --> materiaal matrix Real \in R(VW x VH) met waarden tussen 0 - 1
%q --> uniforme warmte productie per m^2
%Cmet --> coefficient van het metaal (de 1 in de v)
%Cpla --> coefficient van het plastiek (de 0 in de v)   
%BC0,BC1,BC2,BC3 --> bijvoorbeeld
%BC0 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Onder geisoleerde rand
%BC1 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]]; % Rechter 
%BC2 = [['N',1,1,0];['N',2,VB-1,0];['N',VB,VB,0]]; %Boven geisoleerde rand
%BC3 = [['D',1,1,20];['D',2,VH-1,20];['D',VH,VH,20]];% Linker 
% met 'N' neumann, en 'D' diriclet. gevolgd door de 2 indexen van de
% volumes en de vierde en laatste waarde in de rij is waarde van de BC
% zelf.

function [Sol,K] = Harmonic_FVM_NonUniformQ(VW, VH, v, Q, Cmet, Cpla, BC0, BC1, BC2, BC3, p)
H = 1; B = 1; %Hoogte en breedte van het domein
dx = B/(VW-1); dy = H/(VH-1); %Cell grotes

%MateriaalArray
if nargin < 11
    p = 1;
end
MatArray = (1 - v) .^ p * Cpla + v .^ p * Cmet; %verschillende paper
% MatArray = Cmet*simpv + Cpla*(ones(VW,VH)-simpv);

%Aanmaken van matrix en RHS
% K = sparse(VW*VH, VW*VH);
dg   = zeros(VW*VH, 1); % Diagonal
sdg = zeros(VW*VH, 1);  % Subdiagonal
vwdg = zeros(VW*VH, 1); % VB diagonal

%RHS = ones(VW*VH,1);

%temperatuursproductie in de cell
RHS= reshape(Q,[VW*VH,1])*dx*dy;
RHS(1:VW,1) = 1/2*RHS(1:VW,1);%OndersteRij %De border elementen zijn slechts 1/2 of 1/4 de grootte
RHS((VH-1)*VW+1:VH*VW,1) = 1/2*RHS((VH-1)*VW+1:VH*VW,1);%BovensteRij
RHS(1:VW:VW*VH,1) = 1/2*RHS(1:VW:VW*VH,1);%LinkseRij
RHS(VW:VW:VW*VH,1) = 1/2*RHS(VW:VW:VW*VH,1);%RechtseRij

%Stensil voor de inwendige punten 
%Zowel temperatuur en materiaal state zijn cell centered
%Verticale Faces
for i = 1:VW-1 %in breedte
   %Onder
   M = 2/(1/MatArray(i,1)+1/MatArray(i+1,1));%M  = (MatArray(i,1)+MatArray(i+1,1))/2;
   C1 = M*dy/2/dx; C2 = M*dy/2/dx;
   k = i;
   %Vergelijking k
%    K(k,k) = K(k,k) + C1; 
%    K(k,k+1) = K(k,k+1) - C2;
%    %vergelijking k+1
%    K(k+1,k) = K(k+1,k) - C1; 
%    K(k+1,k+1) = K(k+1,k+1) + C2;
   
   dg(k) = dg(k) + C1;
   dg(k+1) = dg(k+1) + C2;
   sdg(k) = sdg(k) - C2;
   
   %Intern
   for j = 2:VH-1 %in hooghte
       M = 2/(1/MatArray(i,j)+1/MatArray(i+1,j));%M = (MatArray(i,j)+MatArray(i+1,j))/2;
       C1 = M*dy/dx; 
       C2 = M*dy/dx;
       k  = i+VW*(j-1);
%        %Vergelijking k
%        K(k,k) = K(k,k) + C1;
%        K(k,k+1) = K(k,k+1) - C2;
%        %vergelijking k+1
%        K(k+1,k) = K(k+1,k) - C1;
%        K(k+1,k+1) = K(k+1,k+1) + C2;
       
       dg(k) = dg(k) + C1;
       dg(k+1) = dg(k+1) + C2;
       sdg(k) = sdg(k) - C2;
   end
   
   %Boven
   M = 2/(1/MatArray(i,VH)+1/MatArray(i+1,VH));%M = (MatArray(i,VH)+MatArray(i+1,VH))/2;
   C1 = M*dy/2/dx; C2 = M*dy/2/dx;
   k = i+VW*(VH-1);
   %Vergelijking k
%    K(k,k) = K(k,k) + C1; K(k,k+1) = K(k,k+1) - C2;
%    %vergelijking k+1
%    K(k+1,k) = K(k+1,k) - C1; K(k+1,k+1) = K(k+1,k+1) + C2;
   
   dg(k) = dg(k) + C1;
   dg(k+1) = dg(k+1) + C2;
   sdg(k) = sdg(k) - C2;
end

%Horizontale Faces
for j = 1:VH-1 %in hooghte
    %Links
    M = 2/(1/MatArray(1,j)+1/MatArray(1,j+1));%M = (MatArray(1,j)+MatArray(1,j+1))/2;
    C1 = M*dx/2/dy; C2 = M*dx/2/dy;
    k = 1+VW*(j-1);
    %Vergelijking k
%     K(k,k) = K(k,k) + C1; K(k,k+VW) = K(k,k+VW) - C2;
%     %vergelijking k+VB
%     K(k+VW,k) = K(k+VW,k) - C1; K(k+VW,k+VW) = K(k+VW,k+VW) + C2;
    
    dg(k) = dg(k) + C1;
    dg(k+VW) = dg(k+VW) + C2;
    vwdg(k) = vwdg(k) - C2;
    %Intern
    for i = 2:VW-1 %in breedte
       M = 2/(1/MatArray(i,j)+1/MatArray(i,j+1));%M = (MatArray(i,j)+MatArray(i,j+1))/2;
       C1 = M*dx/dy; C2 = M*dx/dy;
       k = i+VW*(j-1);
       %Vergelijking k
%        K(k,k) = K(k,k) + C1; K(k,k+VW) = K(k,k+VW) - C2;
%        %vergelijking k+VB
%        K(k+VW,k) = K(k+VW,k) - C1; K(k+VW,k+VW) = K(k+VW,k+VW) + C2;
       
        dg(k) = dg(k) + C1;
        dg(k+VW) = dg(k+VW) + C2;
        vwdg(k) = vwdg(k) - C2;
    end
   
    %Rechts
    M = 2/(1/MatArray(VW,j)+1/MatArray(VW,j+1)); %M = (MatArray(VW,j)+MatArray(VW,j+1))/2;
    C1 = M*dx/2/dy; C2 = M*dx/2/dy;
    k = VW+VW*(j-1);
    %Vergelijking k
%     K(k,k) = K(k,k) + C1; K(k,k+VW) = K(k,k+VW) - C2;
%     %vergelijking k+VB
%     K(k+VW,k) = K(k+VW,k) - C1; K(k+VW,k+VW) = K(k+VW,k+VW) + C2;
    
    dg(k) = dg(k) + C1;
    dg(k+VW) = dg(k+VW) + C2;
    vwdg(k) = vwdg(k) - C2;
end

%spy(K);

%%%%%%%%%%%%%%%BC
%%%%Boundary conditions onder
PW= 10^8; %Penaltywaarde
%Kleine Cell
if BC0(1,1)=='D'
       k = 1;
%        K(k,k) = K(k,k)+ PW;
       RHS(k) = RHS(k) + BC0(1,4)*PW;
       dg(k) = dg(k) + PW;
else %Neumann   
    C1 =  BC0(1,4)*dx/2;
    k = 1;
    RHS(k) = RHS(k) + C1;
end
%%Interne Cellen
for i = 2:size(BC0,1)-1 
    if BC0(i,1)=='D'
        for j = BC0(i,2):BC0(i,3)
            k = j;
%             K(k,k) = K(k,k)+ PW;
            RHS(k) = RHS(k) + BC0(i,4)*PW;
            dg(k) = dg(k) + PW;
        end
    else %Neumann   
        for j = BC0(i,2):BC0(i,3)
        C1 =  BC0(i,4)*dx;
        k = j;
        RHS(k) = RHS(k) + C1;
        end
    end
end
%Eindcell
if BC0(1,1)=='D' %diriclet
       k = VW;
%        K(k,k) = K(k,k)+ PW;
       RHS(k) = RHS(k) + BC0(size(BC0,1),4)*PW;
       dg(k) = dg(k) + PW;
else %Neumann   
    C1 =  BC0(size(BC0,1),4)*dx/2;
    k = VW;
    RHS(k) = RHS(k) + C1;
end

%%%Boundary conditions Rechts %BC1
%Kleine Cell
if BC1(1,1)=='D'
       k = VW;
%        K(k,k) = K(k,k)+ PW;
       RHS(k) = RHS(k) + BC1(1,4)*PW;
       dg(k) = dg(k) + PW;
else %Neumann   
    C1 =  BC1(1,4)*dx/2;
    k = VW;
    RHS(k) = RHS(k) + C1;
end
%%Interne Cellen
for i = 2:size(BC1,1)-1 
    if BC1(i,1)=='D'
        for j = BC1(i,2):BC1(i,3)
            k = j*VW;
%             K(k,k) = K(k,k)+ PW;
            RHS(k) = RHS(k) + BC1(i,4)*PW;
            dg(k) = dg(k) + PW;
        end
    else %Neumann   
        for j = BC1(i,2):BC1(i,3)
            C1 =  BC1(i,4)*dx;
            k = j*VW;
            RHS(k) = RHS(k) + C1;
        end
    end
end
%Eindcell
if BC1(1,1)=='D' %diriclet
       k = VW*VH;
%        K(k,k) = K(k,k)+ PW;
       RHS(k) = RHS(k) + BC1(size(BC1,1),4)*PW;
       dg(k) = dg(k) + PW;
else %Neumann   
    C1 =  BC1(size(BC1,1),4)*dx/2;
    k = VW*VH;
    RHS(k) = RHS(k) + C1;
end

%%%Boundary conditions Boven %BC2
%Kleine Cell
if BC2(1,1)=='D'
       k = 1+(VH-1)*VW;
%        K(k,k) = K(k,k)+ PW;
       RHS(k) = RHS(k) + BC2(1,4)*PW;
       dg(k) = dg(k)+ PW;
else %Neumann   
    C1 =  BC2(1,4)*dx/2;
    k = 1+(VH-1)*VW;
    RHS(k) = RHS(k) + C1;
end
%%Interne Cellen
for i = 2:size(BC2,1)-1 
    if BC2(i,1)=='D'
        for j = BC2(i,2):BC2(i,3)
            k = j+(VH-1)*VW;
%             K(k,k) = K(k,k)+ PW;
            RHS(k) = RHS(k) + BC2(i,4)*PW;
            dg(k) = dg(k) + PW;
        end
    else %Neumann   
        for j = BC2(i,2):BC2(i,3)
            C1 =  BC2(i,4)*dx;
            k = j+(VH-1)*VW;
            RHS(k) = RHS(k) + C1;
        end
    end
end
%Eindcell
if BC2(1,1)=='D' %diriclet
       k = VW*VH;
%        K(k,k) = K(k,k)+ PW;
       RHS(k) = RHS(k) + BC2(size(BC2,1),4)*PW;
       dg(k) = dg(k) + PW;
else %Neumann   
    C1 =  BC2(size(BC2,1),4)*dx/2;
    k = VW*VH;
    RHS(k) = RHS(k) + C1;
end

%%%Boundary conditions Links %BC3
%Kleine Cell
if BC3(1,1)=='D'
       k = 1;
%        K(k,k) = K(k,k)+ PW;
       RHS(k) = RHS(k) + BC3(1,4)*PW;
       dg(k) = dg(k) + PW;
else %Neumann   
    C1 =  BC3(1,4)*dx/2;
    k = 1+(VH-1)*VW;
    RHS(k) = RHS(k) + C1;
end
%%Interne Cellen
for i = 2:size(BC3,1)-1 
    if BC3(i,1)=='D'
        for j = BC3(i,2):BC3(i,3)
            k = 1 + (j-1)*VW;
%             K(k,k) = K(k,k)+ PW;
            dg(k) = dg(k) + PW;
            RHS(k) = RHS(k) + BC3(i,4)*PW;
        end
    else %Neumann   
        for j = BC3(i,2):BC3(i,3)
            C1 =  BC3(i,4)*dx;
            k = 1 + (j-1)*VW;
            RHS(k) = RHS(k) + C1;
        end
    end
end
%Eindcell
if BC3(1,1)=='D' %diriclet
       k = 1+VW*(VH-1);
%        K(k,k) = K(k,k)+ PW;
       RHS(k) = RHS(k) + BC3(size(BC3,1),4)*PW;
       dg(k) = dg(k)+ PW;
else %Neumann   
    C1 =  BC3(size(BC3,1),4)*dx/2;
    k = 1+VW*(VH-1);
    RHS(k) = RHS(k) + C1;
end

%Oplossen Systeem
K = spdiags([vwdg, sdg, dg, [0; sdg(1:end-1)], [zeros(VW, 1); vwdg(1:end-VW)]], [-VW, -1, 0, 1, VW], VW*VH, VW*VH);
Sol = K\RHS;
% Sol = bicgstab(K, RHS);
end