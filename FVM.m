function [Sol,K] = FVM(VW, VH, v, q, Cmet, Cpla)
H = 1; B = 1; %Hoogte en breedte van het domein
dx = B/(VW-1); dy = H/(VH-1); %Cell grotes

%MateriaalArray
MatArray = Cmet*v + Cpla*(ones(VW,VH)-v);

%Aanmaken van matrix en RHS
K = sparse(VW*VH, VW*VH);
RHS = ones(VW*VH,1);

%temperatuursproductie in de cell
RHS= q*dx*dy*RHS;
RHS(1:VW,1) = 1/2*RHS(1:VW,1);%OndersteRij %De border elementen zijn slechts 1/2 of 1/4 de grootte
RHS((VH-1)*VW+1:VH*VW,1) = 1/2*RHS((VH-1)*VW+1:VH*VW,1);%BovensteRij
RHS(1:VW:VW*VH,1) = 1/2*RHS(1:VW:VW*VH,1);%LinkseRij
RHS(VW:VW:VW*VH,1) = 1/2*RHS(VW:VW:VW*VH,1);%RechtseRij

%Stensil voor de inwendige punten 
%Zowel temperatuur en materiaal state zijn cell centered
%Verticale Faces
for i = 1:VW-1 %in breedte
   %Onder
   M = (MatArray(i,1)+MatArray(i+1,1))/2;
   C1 = M*dy/2/dx; C2 = M*dy/2/dx;
   k = i;
   %Vergelijking k
   K(k,k) = K(k,k) + C1; K(k,k+1) = K(k,k+1) - C2;
   %vergelijking k+1
   K(k+1,k) = K(k+1,k) - C1; K(k+1,k+1) = K(k+1,k+1) + C2;
   
   %Intern
   for j = 2:VH-1 %in hooghte
       M = (MatArray(i,j)+MatArray(i+1,j))/2;
       C1 = M*dy/dx; 
       C2 =M*dy/dx;
       k = i+VW*(j-1);
       %Vergelijking k
       K(k,k) = K(k,k) + C1;
       K(k,k+1) = K(k,k+1) - C2;
       %vergelijking k+1
       K(k+1,k) = K(k+1,k) - C1;
       K(k+1,k+1) = K(k+1,k+1) + C2;
   end
   
   %Boven
   M = (MatArray(i,VH)+MatArray(i+1,VH))/2;
   C1 = M*dy/2/dx; C2 = M*dy/2/dx;
   k = i+VW*(VH-1);
   %Vergelijking k
   K(k,k) = K(k,k) + C1; K(k,k+1) = K(k,k+1) - C2;
   %vergelijking k+1
   K(k+1,k) = K(k+1,k) - C1; K(k+1,k+1) = K(k+1,k+1) + C2;
end

%Horizontale Faces
for j = 1:VH-1 %in hooghte
    %Links
    M = (MatArray(1,j)+MatArray(1,j+1))/2;
    C1 = M*dx/2/dy; C2 = M*dx/2/dy;
    k = 1+VW*(j-1);
    %Vergelijking k
    K(k,k) = K(k,k) + C1; K(k,k+VW) = K(k,k+VW) - C2;
    %vergelijking k+VB
    K(k+VW,k) = K(k+VW,k) - C1; K(k+VW,k+VW) = K(k+VW,k+VW) + C2;
    
    %Intern
    for i = 2:VW-1 %in breedte
       M = (MatArray(i,j)+MatArray(i,j+1))/2;
       C1 = M*dx/dy; C2 = M*dx/dy;
       k = i+VW*(j-1);
       %Vergelijking k
       K(k,k) = K(k,k) + C1; K(k,k+VW) = K(k,k+VW) - C2;
       %vergelijking k+VB
       K(k+VW,k) = K(k+VW,k) - C1; K(k+VW,k+VW) = K(k+VW,k+VW) + C2;
    end
   
    %Rechts
    M = (MatArray(VW,j)+MatArray(VW,j+1))/2;
    C1 = M*dx/2/dy; C2 = M*dx/2/dy;
    k = VW+VW*(j-1);
    %Vergelijking k
    K(k,k) = K(k,k) + C1; K(k,k+VW) = K(k,k+VW) - C2;
    %vergelijking k+VB
    K(k+VW,k) = K(k+VW,k) - C1; K(k+VW,k+VW) = K(k+VW,k+VW) + C2;
end

%spy(K);

%%%%%%%%%%%%%%%BC
%%%%Neumann -> niets voor homogene boundaryconditions %Richting in het
%%%%domein
%Linkse Rand
% DT1 = 0;
%Kleine Cell
%C1 = DT1*MatArray(1,1)*dy/2;
%k = 1;
%RHS(k) = RHS(k) + C1;
%Gewone Cell
% for i =2:VH-1
%     C1 = DT1*MatArray(1,i)*dy;
%     k = 1+(i-1)*VB;
%     RHS(k) = RHS(k) + C1;
% end
%KleineCell
%C1 = DT1*MatArray(1,VH)*dy/2;
%k = 1+(VH-1)*VB
%RHS(k) = RHS(k) + C1;

%Rechtse rand 
% DT2 = 0;
%KleineCell
%C1 = DT2*MatArray(VB,1)*dy/2;
%      k = VB;
%      RHS(k) = RHS(k) + C1;
%Gewone Grootte
% for i =2:VH-1
%      C1 = DT2*MatArray(VB,i)*dy;
%      k = i*VB;
%      RHS(k) = RHS(k) + C1;
% end
%Kleine Cell
% C1 = DT2*MatArray(VB,VH)*dy/2;
%      k = VH*VB;
%      RHS(k) = RHS(k) + C1;

%Boven rand
DT3 = 0;
%Kleine Cell
C1 = DT3*dx/2;  % TODO 
k = 1+(VH-1)*VW;
RHS(k) = RHS(k) + C1;
%Gewone Grootte
for i =2:VW-1
    C1 = DT3*dx;
    k = i+(VH-1)*VW;
    RHS(k) = RHS(k) + C1;
end
%Kleine Cell
C1 = DT3*dx/2;
k = VH*VW;
RHS(k) = RHS(k) + C1;


%Onderrand
DT4 = 0;
%Kleine Cell
C1 = DT4*dx/2;
k = 1;
RHS(k) = RHS(k) + C1;
%Gewone Grootte
for i =2:VW-1
    C1 = DT4*dx;
    k = i;
    RHS(k) = RHS(k) + C1;
end
%KleineCell
C1 = DT4*dx/2;
k = VW;
RHS(k) = RHS(k) + C1;

%%%%Diriclet
PW= 10^8; %Penaltywaarde
%links
T1 = 20;
for i = 1:VH
    k = 1 + (i-1)*VW;
    %K(k,:) = zeros(1,VB*VH); 
    K(k,k) = K(k,k)+ PW;RHS(k) = RHS(k) + T1*PW;
end
%rechts
T2 = 20;
for i = 1:VH
    k = i*VW;
    %K(k,:) = zeros(1,VB*VH); 
    K(k,k) = K(k,k)+ PW;RHS(k) = RHS(k) + T2*PW;
end

% %onder
% T3 = 0;
% for i = 1:VB
%     k = i;
%     K(k,:) = zeros(1,VB*VH); K(k,k) = 1;RHS(k) = T3;
% end
% %boven
% T4 = 0;
% for i = 1:VB
%     k = i + VB*(VH-1);
%     K(k,:) = zeros(1,VB*VH); K(k,k) = 1;RHS(k) = T3;
% end

%Oplossen Systeem
Sol = K\RHS;
end