function AG = Adjoint_Gradient(VB,VH,Varray,L,SOL)%moet synchroon blijven met Bovenstaande Methode
H = 1; B = 1; %Hoogte en breedte van het domein
dx = B/(VB-1); dy = H/(VH-1); %Cell grotes
Cmet = 65; %Metaal
Cpla = 0.2; %Plastiek

G = zeros(VB*VH,VB*VH);

%Stensil voor de inwendige punten 
%Zowel temperatuur en materiaal state zijn cell centered
%Verticale Faces

R = sparse(VB*VH,VB*VH);
for i = 1:VB-1 %in breedte
   %Onder
   k = i;
   R1 = dy/dx/4*(Cmet-Cpla);
   R2 = dy/dx/4*(Cmet-Cpla);
   
%    Vergelijking k
%    R(k,k) = R(k,k) + R1; 
%    R(k,k+1) = R(k,k+1) - R2;
%    R(k+1,k) = R(k+1,k) - R1; 
%    R(k+1,k+1) = R(k+1,k+1) + R2;
%    G(:,k) = G(:,k)+R*SOL;
%   
%    R = sparse(VB*VH,VB*VH);
   
   %vergelijking k+1
%    R(k,k) = R(k,k) + R1; 
%    R(k,k+1) = R(k,k+1) - R2;
%    R(k+1,k) = R(k+1,k) - R1;
%    R(k+1,k+1) = R(k+1,k+1) + R2;
%    G(:,k+1) = G(:,k+1)+R*SOL;
   
   G(k,   k) = G(k,   k) + R1 * SOL(k) - R2 * SOL(k+1);
   G(k+1, k) = G(k+1, k) - R1 * SOL(k) + R2 * SOL(k+1);
   
   G(k,   k+1) = G(k,   k+1) + R1 * SOL(k) - R2 * SOL(k+1);
   G(k+1, k+1) = G(k+1, k+1) - R1 * SOL(k) + R2 * SOL(k+1);
   
   R = sparse(VB*VH,VB*VH); 
   
   %Intern
   for j = 2:VH-1 %in hooghte
        k = i+VB*(j-1);
        R1 = dy/dx/2*(Cmet-Cpla); 
        R2 = dy/dx/2*(Cmet-Cpla);
       
%         T = G;
%         
       %Vergelijking k
%        R(k,k) = R(k,k) + R1; 
%        R(k,k+1) = R(k,k+1) - R2;
%        R(k+1,k) = R(k+1,k) - R1; 
%        R(k+1,k+1) = R(k+1,k+1) + R2;
%        G(:,k) = G(:,k)+R*SOL;
%        R = sparse(VB*VH,VB*VH);
%        
%        %vergelijking k+1
%        R(k,k) = R(k,k) + R1;
%        R(k,k+1) = R(k,k+1) - R2;
%        R(k+1,k) = R(k+1,k) - R1;
%        R(k+1,k+1) = R(k+1,k+1) + R2;
%        G(:,k+1) = G(:,k+1)+R*SOL;
%        
%        T(k,   k:k+1) =  R1 * SOL(k) - R2 * SOL(k+1);
%        T(k+1, k:k+1) = -R1 * SOL(k) + R2 * SOL(k+1);
       
       G(k,   k) = G(k,   k) + R1 * SOL(k) - R2 * SOL(k+1);
       G(k+1, k) = G(k+1, k) - R1 * SOL(k) + R2 * SOL(k+1);
   
       G(k,   k+1) = G(k,   k+1) + R1 * SOL(k) - R2 * SOL(k+1);
       G(k+1, k+1) = G(k+1, k+1) - R1 * SOL(k) + R2 * SOL(k+1);
       
%        R = sparse(VB*VH,VB*VH); 
       
   end
   
   %Boven
   k = i+VB*(VH-1);
   R1 = dy/dx/4*(Cmet-Cpla); R2 = dy/dx/4*(Cmet-Cpla);

%    %Vergelijking k
%    R(k,k) = R(k,k) + R1; R(k,k+1) = R(k,k+1) - R2;
%    R(k+1,k) = R(k+1,k) - R1; R(k+1,k+1) = R(k+1,k+1) + R2;
%    G(:,k) = G(:,k)+R*SOL;
%    R = sparse(VB*VH,VB*VH);
% 
%    %vergelijking k+1
%    R(k,k) = R(k,k) + R1; R(k,k+1) = R(k,k+1) - R2;
%    R(k+1,k) = R(k+1,k) - R1; R(k+1,k+1) = R(k+1,k+1) + R2;
%    G(:,k+1) = G(:,k+1)+R*SOL;
%    R = sparse(VB*VH,VB*VH);
   
    G(k,   k) = G(k,   k) + R1 * SOL(k) - R2 * SOL(k+1);
    G(k+1, k) = G(k+1, k) - R1 * SOL(k) + R2 * SOL(k+1);
   
    G(k,   k+1) = G(k,   k+1) + R1 * SOL(k) - R2 * SOL(k+1);
    G(k+1, k+1) = G(k+1, k+1) - R1 * SOL(k) + R2 * SOL(k+1);
end

%Horizontale Faces
for j = 1:VH-1 %in hooghte
    %Links
    k = 1+VB*(j-1);
    R1 = dx/dy/4*(Cmet-Cpla); R2 = dx/dy/4*(Cmet-Cpla);
    
    %Vergelijking k
%     R(k,k) = R(k,k) + R1; 
%     R(k,k+VB) = R(k,k+VB) - R2;
%     R(k+VB,k) = R(k+VB,k) - R1; 
%     R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
%     
%     G(:,k) = G(:,k)+R*SOL;
%     R = sparse(VB*VH,VB*VH);
%     
%     %vergelijking k+VB
%     R(k,k) = R(k,k) + R1; R(k,k+VB) = R(k,k+VB) - R2;
%     R(k+VB,k) = R(k+VB,k) - R1; R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
%     G(:,k+VB) = G(:,k+VB)+R*SOL;
%     R = sparse(VB*VH,VB*VH);
    
    G(k,    k) = G(k,    k) + R1 * SOL(k) - R2 * SOL(k+VB);
    G(k+VB, k) = G(k+VB, k) - R1 * SOL(k) + R2 * SOL(k+VB);
   
    G(k,    k+VB) = G(k,    k+VB) + R1 * SOL(k) - R2 * SOL(k+VB);
    G(k+VB, k+VB) = G(k+VB, k+VB) - R1 * SOL(k) + R2 * SOL(k+VB);
    
    %Intern
    for i = 2:VB-1 %in breedte
        k = i+VB*(j-1);
        R1 = dx/dy/2*(Cmet-Cpla); R2 = dx/dy/2*(Cmet-Cpla);
%         
%         %Vergelijking k
%         R(k,k) = R(k,k) + R1; R(k,k+VB) = R(k,k+VB) - R2;
%         R(k+VB,k) = R(k+VB,k) - R1; R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
%         G(:,k) = G(:,k)+R*SOL;
%         R = sparse(VB*VH,VB*VH);
%         
%         %vergelijking k+VB
%         R(k,k) = R(k,k) + R1; R(k,k+VB) = R(k,k+VB) - R2;
%         R(k+VB,k) = R(k+VB,k) - R1; R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
%         G(:,k+VB) = G(:,k+VB)+R*SOL;
%         R = sparse(VB*VH,VB*VH);
        
        G(k,    k) = G(k,    k) + R1 * SOL(k) - R2 * SOL(k+VB);
        G(k+VB, k) = G(k+VB, k) - R1 * SOL(k) + R2 * SOL(k+VB);
   
        G(k,    k+VB) = G(k,    k+VB) + R1 * SOL(k) - R2 * SOL(k+VB);
        G(k+VB, k+VB) = G(k+VB, k+VB) - R1 * SOL(k) + R2 * SOL(k+VB);
    end
   
    %Rechts
    k = VB+VB*(j-1);
    R1 = dx/dy/4*(Cmet-Cpla); R2 = dx/dy/4*(Cmet-Cpla);
% 
%     %Vergelijking k
%     R(k,k) = R(k,k)+ R1; 
%     R(k,k+VB) = R(k,k+VB) - R2;
%     R(k+VB,k) = R(k+VB,k) - R1; 
%     R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
%     G(:,k) = G(:,k)+R*SOL;
%     R = sparse(VB*VH,VB*VH);
%     
%     %vergelijking k+VB
%     R(k,k) = R(k,k) + R1; 
%     R(k,k+VB) = R(k,k+VB) - R2;
%     R(k+VB,k) = R(k+VB,k) - R1; 
%     R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
%     G(:,k+VB) = G(:,k+VB)+R*SOL;
%     R = sparse(VB*VH,VB*VH);

      G(k,    k) = G(k,    k) + R1 * SOL(k) - R2 * SOL(k+VB);
      G(k+VB, k) = G(k+VB, k) - R1 * SOL(k) + R2 * SOL(k+VB);
   
      G(k,    k+VB) = G(k,    k+VB) + R1 * SOL(k) - R2 * SOL(k+VB);
      G(k+VB, k+VB) = G(k+VB, k+VB) - R1 * SOL(k) + R2 * SOL(k+VB);
end

AG = L'*G;
end
