function AG = Adjoint_Gradient(VB,VH,Varray,L,SOL)%moet synchroon blijven met Bovenstaande Methode
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

   %ADJ
   R1 = dy/dx/4*(Cmet-Cpla);
   R2 = dy/dx/4*(Cmet-Cpla);
   %Vergelijking k
   R(k,k) = R(k,k) + R1; 
   R(k,k+1) = R(k,k+1) - R2;
   R(k+1,k) = R(k+1,k) - R1; 
   R(k+1,k+1) = R(k+1,k+1) + R2;
   G(:,k) = G(:,k)+R*SOL;
   R = sparse(VB*VH,VB*VH);
   R(k,k) = R(k,k) + R1;
   R(k,k+1) = R(k,k+1) - R2;
   %vergelijking k+1
   
   R(k+1,k) = R(k+1,k) - R1; 
   R(k+1,k+1) = R(k+1,k+1) + R2;
   G(:,k+1) = G(:,k+1)+R*SOL;
   R = sparse(VB*VH,VB*VH); 
   
   %Intern
   for j = 2:VH-1 %in hooghte
        k = i+VB*(j-1);
       %ADJ
       R1 = dy/dx/2*(Cmet-Cpla); R2 = dy/dx/2*(Cmet-Cpla);
       %Vergelijking k
       R(k,k) = R(k,k) + R1;
       R(k,k+1) = R(k,k+1) - R2;
       %vergelijking k+1
       R(k+1,k) = R(k+1,k) - R1;
       R(k+1,k+1) = R(k+1,k+1) + R2;
       G(:,k) = G(:,k)+R*SOL;
       R = sparse(VB*VH,VB*VH);
       
       R(k,k) = R(k,k) + R1;R(k,k+1) = R(k,k+1) - R2;
       R(k+1,k) = R(k+1,k) - R1; R(k+1,k+1) = R(k+1,k+1) + R2;
       G(:,k+1) = G(:,k+1)+R*SOL;
       R = sparse(VB*VH,VB*VH); 
       
   end
   
   %Boven
    k = i+VB*(VH-1);

%ADJ
       R1 = dy/dx/4*(Cmet-Cpla); R2 = dy/dx/4*(Cmet-Cpla);
       %Vergelijking k
       R(k,k) = R(k,k) + R1;
       R(k,k+1) = R(k,k+1) - R2;
%               R(k,k,k) = R(k,k,k) + R1; R(k,k+1,k) = R(k,k+1,k) - R2;
       %vergelijking k+1
       R(k+1,k) = R(k+1,k) - R1;
       R(k+1,k+1) = R(k+1,k+1) + R2;
       G(:,k) = G(:,k)+R*SOL;
       R = sparse(VB*VH,VB*VH);
       R(k,k) = R(k,k) + R1; R(k,k+1) = R(k,k+1) - R2;
%               R(k,k,k+1) = R(k,k,k+1) + R1; R(k,k+1,k+1) = R(k,k+1,k+1) - R2;

       R(k+1,k) = R(k+1,k) - R1; R(k+1,k+1) = R(k+1,k+1) + R2;
       G(:,k+1) = G(:,k+1)+R*SOL;
       R = sparse(VB*VH,VB*VH);
       

       %vergelijking k+1
%        R(k+1,k,k) = R(k+1,k,k) - R1; R(k+1,k+1,k) = R(k+1,k+1,k) + R2;
%        R(k+1,k,k+1) = R(k+1,k,k+1) - R1; R(k+1,k+1,k+1) = R(k+1,k+1,k+1) + R2;
end

%Horizontale Faces
for j = 1:VH-1 %in hooghte
    %Links
    k = 1+VB*(j-1);
    R1 = dx/dy/4*(Cmet-Cpla); R2 = dx/dy/4*(Cmet-Cpla);
    %Vergelijking k
    R(k,k) = R(k,k) + R1;
    R(k,k+VB) = R(k,k+VB) - R2;
    
    %vergelijking k+VB
    R(k+VB,k) = R(k+VB,k) - R1;
    R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
    G(:,k) = G(:,k)+R*SOL;
    R = sparse(VB*VH,VB*VH);
    
    R(k,k) = R(k,k) + R1; R(k,k+VB) = R(k,k+VB) - R2;
    R(k+VB,k) = R(k+VB,k) - R1; R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
    G(:,k+VB) = G(:,k+VB)+R*SOL;
    R = sparse(VB*VH,VB*VH);
    %Intern
    for i = 2:VB-1 %in breedte
       k = i+VB*(j-1);
        R1 = dx/dy/2*(Cmet-Cpla); R2 = dx/dy/2*(Cmet-Cpla);
        %Vergelijking k
    R(k,k) = R(k,k) + R1;
    R(k,k+VB) = R(k,k+VB) - R2;

    %vergelijking k+VB
    R(k+VB,k) = R(k+VB,k) - R1;
    R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
    G(:,k) = G(:,k)+R*SOL;
    R = sparse(VB*VH,VB*VH);
    
    R(k,k) = R(k,k) + R1; R(k,k+VB) = R(k,k+VB) - R2;
    R(k+VB,k) = R(k+VB,k) - R1; R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
    G(:,k+VB) = G(:,k+VB)+R*SOL;
    R = sparse(VB*VH,VB*VH);
%     R(k,k,k) = R(k,k,k) + R1; R(k,k+VB,k) = R(k,k+VB,k) - R2;
%     R(k+VB,k,k) = R(k+VB,k,k) - R1; R(k+VB,k+VB,k) = R(k+VB,k+VB,k) + R2;

%     R(k,k,k+VB) = R(k,k,k+VB) + R1; R(k,k+VB,k+VB) = R(k,k+VB,k+VB) - R2;
%     R(k+VB,k,k+VB) = R(k+VB,k,k+VB) - R1; R(k+VB,k+VB,k+VB) = R(k+VB,k+VB,k+VB) + R2;

    end
   
    %Rechts
    k = VB+VB*(j-1);
    R1 = dx/dy/4*(Cmet-Cpla); R2 = dx/dy/4*(Cmet-Cpla);
    %Vergelijking k
    R(k,k) = R(k,k) + R1;
    R(k,k+VB) = R(k,k+VB) - R2;
    
    %vergelijking k+VB
    R(k+VB,k) = R(k+VB,k) - R1;
    R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
    G(:,k) = G(:,k)+R*SOL;
    R = sparse(VB*VH,VB*VH);
    
    R(k,k) = R(k,k) + R1; R(k,k+VB) = R(k,k+VB) - R2;
    R(k+VB,k) = R(k+VB,k) - R1; R(k+VB,k+VB) = R(k+VB,k+VB) + R2;
    G(:,k+VB) = G(:,k+VB)+R*SOL;
    R = sparse(VB*VH,VB*VH);
end

AG = L'*G;
end