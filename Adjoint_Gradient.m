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
K = zeros(VB*VH,VB*VH);
G = zeros(VB*VH,VB*VH);

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
    R(k,k,k+VB) = R(k,k,k+VB) + R1; K(k,k+VB,k+VB) = R(k,k+VB,k+VB) - R2;
    %vergelijking k+VB
    R(k+VB,k,k) = R(k+VB,k,k) - R1; R(k+VB,k+VB,k) = R(k+VB,k+VB,k) + R2;
    R(k+VB,k,k+VB) = R(k+VB,k,k+VB) - R1; R(k+VB,k+VB,k+VB) = R(k+VB,k+VB,k+VB) + R2;
end

for i = 1:VH*VB
    G(:,i) = G(:,i) + R(:,:,i)*SOL;
end

AG = L'*G;
end