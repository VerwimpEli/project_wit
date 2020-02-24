clear;

%Constanten %Mesh wordt volledig vierkant verondersteld %Speciale
%Vierkanten aan de randen zodat de berekende temperaturen het volledige
%domein insluiten
VH = 30; VB = 30; % Aantal volumes in de hoogte en breedte. Incluisief de kleinere op de randen
rng(500); %Reproduceerbaarheid
Varray = rand(VB,VH);
%Varray = 0.5*ones(VB,VH);
[Sol,K] = FVM(VB,VH,Varray);
% Best niet veel hoger dan 125 --> of spaarse matrices gebruiken

%Visualisatie %Heel ruw komt niet direct overeen met echte systeem
SOL = reshape(Sol,[VB,VH]);
figure(1);
surf(SOL); 
xlabel("X"); ylabel("Y"); zlabel("Temperatuur")

%Benaderen Jacobiaan
% J = FD_J(VB,VH,Varray);
%G_FD = ones(1,VB*VH)*J';
% G_FDT = ones(1,VB*VH)*J;%Of moet dit J transpose zijn?

%Adjoint
L = (K')\-ones(VB*VH,1);
% AG = Adjoint_Gradient(VB,VH,Varray,L,Sol);
AG2 = Adjoint_Gradient2(VB,VH,Varray,L,Sol);

%norm(AG-G_FD) %wss fout
% norm(AG-AG2)

% ERR1 = reshape(AG-AG2,[VB,VH]);


figure(2);
% subplot(1,3,1); 
surf(reshape(AG2,[VB,VH])); 
% title("AG2");
% subplot(1,3,2)
% surf(ERR1); 
% title("ERR1");
% xlabel("X"); ylabel("Y"); zlabel("ERROR1")
% subplot(1,3,3)
% surf(reshape(AG,[VB,VH])); 
% title("AG");


%%%Bereken van de gradient DMV Finite difference
function J = FD_J(VB,VH,Varray)
    Delta = 10^-6;
    Size = size(Varray,1)*size(Varray,2);
    J = zeros(Size,Size);
    Sol = FVM(VB,VH,Varray);
    
    %Loop
    for i = 1:size(Varray,1)
        for j = 1:size(Varray,2)
            Varray2 = Varray; Varray2(i,j) = Varray2(i,j)*(1+Delta);
            FD_Sol = FVM(VB,VH,Varray2);
            J(:,i+ (j-1)*VB )= (FD_Sol - Sol)/(Delta*Varray(i,j)); %of is dit de gradient?
        end
    end
end

function [Sol,K] = FVM(VB,VH,Varray)
H = 1; B = 1; %Hoogte en breedte van het domein
dx = B/(VB-1); dy = H/(VH-1); %Cell grotes
q = 200;
Cmet = 65; %Metaal
Cpla = 0.2; %Plastiek

%MateriaalArray
MatArray = Cmet*Varray + Cpla*(ones(VB,VH)-Varray);

%Aanmaken van matrix en RHS
K = zeros(VB*VH);
RHS = ones(VB*VH,1);

%temperatuursproductie in de cell
RHS= q*dx*dy*RHS;
RHS(1:VB,1) = 1/2*RHS(1:VB,1);%OndersteRij %De border elementen zijn slecht 1/2 of 1/4 de grootte
RHS((VH-1)*VB+1:VH*VB,1) = 1/2*RHS((VH-1)*VB+1:VH*VB,1);%BovensteRij
RHS(1:VB:VB*VH,1) = 1/2*RHS(1:VB:VB*VH,1);%LinkseRij
RHS(VB:VB:VB*VH,1) = 1/2*RHS(VB:VB:VB*VH,1);%RechtseRij

%Stensil voor de inwendige punten 
%Zowel temperatuur en materiaal state zijn cell centered
%Verticale Faces
for i = 1:VB-1 %in breedte
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
       C1 = M*dy/dx; C2 =M*dy/dx;
       k = i+VB*(j-1);
       %Vergelijking k
       K(k,k) = K(k,k) + C1; K(k,k+1) = K(k,k+1) - C2;
       %vergelijking k+1
       K(k+1,k) = K(k+1,k) - C1; K(k+1,k+1) = K(k+1,k+1) + C2;
   end
   
   %Boven
   M = (MatArray(i,VH)+MatArray(i+1,VH))/2;
   C1 = M*dy/2/dx; C2 = M*dy/2/dx;
   k = i+VB*(VH-1);
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
    k = 1+VB*(j-1);
    %Vergelijking k
    K(k,k) = K(k,k) + C1; K(k,k+VB) = K(k,k+VB) - C2;
    %vergelijking k+VB
    K(k+VB,k) = K(k+VB,k) - C1; K(k+VB,k+VB) = K(k+VB,k+VB) + C2;
    
    %Intern
    for i = 2:VB-1 %in breedte
       M = (MatArray(i,j)+MatArray(i,j+1))/2;
       C1 = M*dx/dy; C2 = M*dx/dy;
       k = i+VB*(j-1);
       %Vergelijking k
       K(k,k) = K(k,k) + C1; K(k,k+VB) = K(k,k+VB) - C2;
       %vergelijking k+VB
       K(k+VB,k) = K(k+VB,k) - C1; K(k+VB,k+VB) = K(k+VB,k+VB) + C2;
    end
   
    %Rechts
    M = (MatArray(VB,j)+MatArray(VB,j+1))/2;
    C1 = M*dx/2/dy; C2 = M*dx/2/dy;
    k = VB+VB*(j-1);
    %Vergelijking k
    K(k,k) = K(k,k) + C1; K(k,k+VB) = K(k,k+VB) - C2;
    %vergelijking k+VB
    K(k+VB,k) = K(k+VB,k) - C1; K(k+VB,k+VB) = K(k+VB,k+VB) + C2;
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
C1 = DT3*dx/2;
k = 1+(VH-1)*VB;
RHS(k) = RHS(k) + C1;
%Gewone Grootte
for i =2:VB-1
    C1 = DT3*dx;
    k = i+(VH-1)*VB;
    RHS(k) = RHS(k) + C1;
end
%Kleine Cell
C1 = DT3*dx/2;
k = VH*VB;
RHS(k) = RHS(k) + C1;


%Onderrand
DT4 = 0;
%Kleine Cell
C1 = DT4*dx/2;
k = 1;
RHS(k) = RHS(k) + C1;
%Gewone Grootte
for i =2:VB-1
    C1 = DT4*dx;
    k = i;
    RHS(k) = RHS(k) + C1;
end
%KleineCell
C1 = DT4*dx/2;
k = VB;
RHS(k) = RHS(k) + C1;

%%%%Diriclet
PW= 10^8; %Penaltywaarde
%links
T1 = 0;
for i = 1:VH
    k = 1 + (i-1)*VB;
    %K(k,:) = zeros(1,VB*VH); 
    K(k,k) = K(k,k)+ PW;RHS(k) = RHS(k) + T1*PW;
end
%rechts
T2 = 0;
for i = 1:VH
    k = i*VB;
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

function AG = Adjoint_Gradient(VB,VH,Varray,L,SOL)%moet synchroon blijven met Bovenstaande Methode
H = 1; B = 1; %Hoogte en breedte van het domein
dx = B/(VB-1); dy = H/(VH-1); %Cell grotes
Cmet = 65; %Metaal
Cpla = 0.2; %Plastiek

%MateriaalArray
MatArray = Cmet*Varray + Cpla*(ones(VB,VH)-Varray);
%MatArray = Cmet*ones(VB,VH);
%MatArray(VB/2:VB,:) = MatArray(VB/2:VB,:)/Cmet*Cpla;

%Rbeta
R = zeros(VB*VH,VB*VH,VB*VH);%Zeer Spaarse Tensor
% R = sparse(VB*VH,VB*VH);
K = sparse(VB*VH,VB*VH);
G = sparse(VB*VH,VB*VH);

%Stensil voor de inwendige punten 
%Zowel temperatuur en materiaal state zijn cell centered
%Verticale Faces
for i = 1:VB-1 %in breedte
   %Onder
%    M = (MatArray(i,1)+MatArray(i+1,1))/2;
%    C1 = M*dy/2/dx; C2 = M*dy/2/dx;
    k = i;
%    %Vergelijking k
%    K(k,k) = K(k,k) + C1; K(k,k+1) = K(k,k+1) - C2;
%    %vergelijking k+1
%    K(k+1,k) = K(k+1,k) - C1; K(k+1,k+1) = K(k+1,k+1) + C2;
   
   %ADJ
   R1 = dy/dx/4*(Cmet-Cpla); R2 = dy/dx/4*(Cmet-Cpla);
   %Vergelijking k
   R(k,k,k) = R(k,k,k) + R1; R(k,k+1,k) = R(k,k+1,k) - R2;
   R(k,k,k+1) = R(k,k,k+1) + R1; R(k,k+1,k+1) = R(k,k+1,k+1) - R2;
   %vergelijking k+1
   R(k+1,k,k) = R(k+1,k,k) - R1; R(k+1,k+1,k) = R(k+1,k+1,k) + R2;
   R(k+1,k,k+1) = R(k+1,k,k+1) - R1; R(k+1,k+1,k+1) = R(k+1,k+1,k+1) + R2;
   
   %Intern
   for j = 2:VH-1 %in hooghte
%        M = (MatArray(i,j)+MatArray(i+1,j))/2;
%        C1 = M*dy/dx; C2 =M*dy/dx;
        k = i+VB*(j-1);
%        %Vergelijking k
%        K(k,k) = K(k,k) + C1; K(k,k+1) = K(k,k+1) - C2;
%        %vergelijking k+1
%        K(k+1,k) = K(k+1,k) - C1; K(k+1,k+1) = K(k+1,k+1) + C2;
       
       %ADJ
       R1 = dy/dx/2*(Cmet-Cpla); R2 = dy/dx/2*(Cmet-Cpla);
       %Vergelijking k
       R(k,k,k) = R(k,k,k) + R1; R(k,k+1,k) = R(k,k+1,k) - R2;
       R(k,k,k+1) = R(k,k,k+1) + R1; R(k,k+1,k+1) = R(k,k+1,k+1) - R2;
       %vergelijking k+1
       R(k+1,k,k) = R(k+1,k,k) - R1; R(k+1,k+1,k) = R(k+1,k+1,k) + R2;
       R(k+1,k,k+1) = R(k+1,k,k+1) - R1; R(k+1,k+1,k+1) = R(k+1,k+1,k+1) + R2;
   end
   
   %Boven
%    M = (MatArray(i,VH)+MatArray(i+1,VH))/2;
%    C1 = M*dy/2/dx; C2 = M*dy/2/dx;
    k = i+VB*(VH-1);
%    %Vergelijking k
%    K(k,k) = K(k,k) + C1; K(k,k+1) = K(k,k+1) - C2;
%    %vergelijking k+1
%    K(k+1,k) = K(k+1,k) - C1; K(k+1,k+1) = K(k+1,k+1) + C2;

%ADJ
       R1 = dy/dx/4*(Cmet-Cpla); R2 = dy/dx/4*(Cmet-Cpla);
       %Vergelijking k
       R(k,k,k) = R(k,k,k) + R1; R(k,k+1,k) = R(k,k+1,k) - R2;
       R(k,k,k+1) = R(k,k,k+1) + R1; R(k,k+1,k+1) = R(k,k+1,k+1) - R2;
       %vergelijking k+1
       R(k+1,k,k) = R(k+1,k,k) - R1; R(k+1,k+1,k) = R(k+1,k+1,k) + R2;
       R(k+1,k,k+1) = R(k+1,k,k+1) - R1; R(k+1,k+1,k+1) = R(k+1,k+1,k+1) + R2;
end

%Horizontale Faces
for j = 1:VH-1 %in hooghte
    %Links
%     M = (MatArray(1,j)+MatArray(1,j+1))/2;
%     C1 = M*dx/2/dy; C2 = M*dx/2/dy;
    k = 1+VB*(j-1);
    %Vergelijking k
%     K(k,k) = K(k,k) + C1; K(k,k+VB) = K(k,k+VB) - C2;
%     %vergelijking k+VB
%     K(k+VB,k) = K(k+VB,k) - C1; K(k+VB,k+VB) = K(k+VB,k+VB) + C2;
    R1 = dx/dy/4*(Cmet-Cpla); R2 = dx/dy/4*(Cmet-Cpla);
    %Vergelijking k
    R(k,k,k) = R(k,k,k) + R1; R(k,k+VB,k) = R(k,k+VB,k) - R2;
    R(k,k,k+VB) = R(k,k,k+VB) + R1; R(k,k+VB,k+VB) = R(k,k+VB,k+VB) - R2;
    %vergelijking k+VB
    R(k+VB,k,k) = R(k+VB,k,k) - R1; R(k+VB,k+VB,k) = R(k+VB,k+VB,k) + R2;
    R(k+VB,k,k+VB) = R(k+VB,k,k+VB) - R1; R(k+VB,k+VB,k+VB) = R(k+VB,k+VB,k+VB) + R2;
    %Intern
    for i = 2:VB-1 %in breedte
%        M = (MatArray(i,j)+MatArray(i,j+1))/2;
%        C1 = M*dx/dy; C2 = M*dx/dy;
       k = i+VB*(j-1);
%        %Vergelijking k
%        K(k,k) = K(k,k) + C1; K(k,k+VB) = K(k,k+VB) - C2;
%        %vergelijking k+VB
%        K(k+VB,k) = K(k+VB,k) - C1; K(k+VB,k+VB) = K(k+VB,k+VB) + C2;
        R1 = dx/dy/2*(Cmet-Cpla); R2 = dx/dy/2*(Cmet-Cpla);
        %Vergelijking k
    R(k,k,k) = R(k,k,k) + R1; R(k,k+VB,k) = R(k,k+VB,k) - R2;
    R(k,k,k+VB) = R(k,k,k+VB) + R1; R(k,k+VB,k+VB) = R(k,k+VB,k+VB) - R2;
    %vergelijking k+VB
    R(k+VB,k,k) = R(k+VB,k,k) - R1; R(k+VB,k+VB,k) = R(k+VB,k+VB,k) + R2;
    R(k+VB,k,k+VB) = R(k+VB,k,k+VB) - R1; R(k+VB,k+VB,k+VB) = R(k+VB,k+VB,k+VB) + R2;
    end
   
    %Rechts
%     M = (MatArray(VB,j)+MatArray(VB,j+1))/2;
%     C1 = M*dx/2/dy; C2 = M*dx/2/dy;
    k = VB+VB*(j-1);
    %Vergelijking k
%     K(k,k) = K(k,k) + C1; K(k,k+VB) = K(k,k+VB) - C2;
%     %vergelijking k+VB
%     K(k+VB,k) = K(k+VB,k) - C1; K(k+VB,k+VB) = K(k+VB,k+VB) + C2;
    R1 = dx/dy/4*(Cmet-Cpla); R2 = dx/dy/4*(Cmet-Cpla);
    %Vergelijking k
    R(k,k,k) = R(k,k,k) + R1; R(k,k+VB,k) = R(k,k+VB,k) - R2;
    R(k,k,k+VB) = R(k,k,k+VB) + R1; R(k,k+VB,k+VB) = R(k,k+VB,k+VB) - R2;
    %vergelijking k+VB
    R(k+VB,k,k) = R(k+VB,k,k) - R1; R(k+VB,k+VB,k) = R(k+VB,k+VB,k) + R2;
    R(k+VB,k,k+VB) = R(k+VB,k,k+VB) - R1; R(k+VB,k+VB,k+VB) = R(k+VB,k+VB,k+VB) + R2;
end

for i = 1:VH*VB
    G(:,i) = G(:,i) + R(:,:,i)*SOL;
end

AG = L'*G;
end

function AG = Adjoint_Gradient2(VB,VH,Varray,L,SOL)%moet synchroon blijven met Bovenstaande Methode
H = 1; B = 1; %Hoogte en breedte van het domein
dx = B/(VB-1); dy = H/(VH-1); %Cell grotes
Cmet = 65; %Metaal
Cpla = 0.2; %Plastiek

G = sparse(VB*VH,VB*VH);

%Stensil voor de inwendige punten 
%Zowel temperatuur en materiaal state zijn cell centered
%Verticale Faces

R = sparse(VB*VH,VB*VH);
for i = 1:VB-1 %in breedte
   %Onder
   k = i;
   R1 = dy/dx/4*(Cmet-Cpla); R2 = dy/dx/4*(Cmet-Cpla);
   
   %Vergelijking k
   R(k,k) = R(k,k) + R1; R(k,k+1) = R(k,k+1) - R2;
   R(k+1,k) = R(k+1,k) - R1; R(k+1,k+1) = R(k+1,k+1) + R2;
   G(:,k) = G(:,k)+R*SOL;
   R = sparse(VB*VH,VB*VH);
   
   %vergelijking k+1
   R(k,k) = R(k,k) + R1; R(k,k+1) = R(k,k+1) - R2;
   R(k+1,k) = R(k+1,k) - R1; R(k+1,k+1) = R(k+1,k+1) + R2;
   G(:,k+1) = G(:,k+1)+R*SOL;
   R = sparse(VB*VH,VB*VH); 
   
   %Intern
   for j = 2:VH-1 %in hooghte
       k = i+VB*(j-1);
       R1 = dy/dx/2*(Cmet-Cpla); R2 = dy/dx/2*(Cmet-Cpla);
       
       %Vergelijking k
       R(k,k) = R(k,k) + R1; R(k,k+1) = R(k,k+1) - R2;
       R(k+1,k) = R(k+1,k) - R1; R(k+1,k+1) = R(k+1,k+1) + R2;
       G(:,k) = G(:,k)+R*SOL;
       R = sparse(VB*VH,VB*VH);
       
       %vergelijking k+1
       R(k,k) = R(k,k) + R1;R(k,k+1) = R(k,k+1) - R2;
       R(k+1,k) = R(k+1,k) - R1; R(k+1,k+1) = R(k+1,k+1) + R2;
       G(:,k+1) = G(:,k+1)+R*SOL;
       R = sparse(VB*VH,VB*VH); 
   end
   
   %Boven
   k = i+VB*(VH-1);
   R1 = dy/dx/4*(Cmet-Cpla); R2 = dy/dx/4*(Cmet-Cpla);

   %Vergelijking k
   R(k,k) = R(k,k) + R1; R(k,k+1) = R(k,k+1) - R2;
   R(k+1,k) = R(k+1,k) - R1; R(k+1,k+1) = R(k+1,k+1) + R2;
   G(:,k) = G(:,k)+R*SOL;
   R = sparse(VB*VH,VB*VH);

   %vergelijking k+1
   R(k,k) = R(k,k) + R1; R(k,k+1) = R(k,k+1) - R2;
   R(k+1,k) = R(k+1,k) - R1; R(k+1,k+1) = R(k+1,k+1) + R2;
   G(:,k+1) = G(:,k+1)+R*SOL;
   R = sparse(VB*VH,VB*VH);
end

%Horizontale Faces
for j = 1:VH-1 %in hooghte
    %Links
    k = 1+VB*(j-1);
    R1 = dx/dy/4*(Cmet-Cpla); R2 = dx/dy/4*(Cmet-Cpla);
    
    %Vergelijking k
    R(k,k) = R(k,k) + R1; R(k,k+VB) = R(k,k+VB) - R2;
    R(k+VB,k) = R(k+VB,k) - R1; R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
    G(:,k) = G(:,k)+R*SOL;
    R = sparse(VB*VH,VB*VH);
    
    %vergelijking k+VB
    R(k,k) = R(k,k) + R1; R(k,k+VB) = R(k,k+VB) - R2;
    R(k+VB,k) = R(k+VB,k) - R1; R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
    G(:,k+VB) = G(:,k+VB)+R*SOL;
    R = sparse(VB*VH,VB*VH);
    
    %Intern
    for i = 2:VB-1 %in breedte
        k = i+VB*(j-1);
        R1 = dx/dy/2*(Cmet-Cpla); R2 = dx/dy/2*(Cmet-Cpla);
        
        %Vergelijking k
        R(k,k) = R(k,k) + R1; R(k,k+VB) = R(k,k+VB) - R2;
        R(k+VB,k) = R(k+VB,k) - R1; R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
        G(:,k) = G(:,k)+R*SOL;
        R = sparse(VB*VH,VB*VH);
        
        %vergelijking k+VB
        R(k,k) = R(k,k) + R1; R(k,k+VB) = R(k,k+VB) - R2;
        R(k+VB,k) = R(k+VB,k) - R1; R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
        G(:,k+VB) = G(:,k+VB)+R*SOL;
        R = sparse(VB*VH,VB*VH);
    end
   
    %Rechts
    k = VB+VB*(j-1);
    R1 = dx/dy/4*(Cmet-Cpla); R2 = dx/dy/4*(Cmet-Cpla);

    %Vergelijking k
    R(k,k) = R(k,k) + R1; R(k,k+VB) = R(k,k+VB) - R2;
    R(k+VB,k) = R(k+VB,k) - R1; R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
    G(:,k) = G(:,k)+R*SOL;
    R = sparse(VB*VH,VB*VH);
    
    %vergelijking k+VB
    R(k,k) = R(k,k) + R1; R(k,k+VB) = R(k,k+VB) - R2;
    R(k+VB,k) = R(k+VB,k) - R1; R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
    G(:,k+VB) = G(:,k+VB)+R*SOL;
    R = sparse(VB*VH,VB*VH);
end

AG = L'*G;
end
