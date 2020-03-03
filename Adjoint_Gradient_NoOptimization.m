function AG = Adjoint_Gradient_NoOptimization(VB,VH,v,L,SOL)%moet synchroon blijven met Bovenstaande Methode
H = 1; B = 1; %Hoogte en breedte van het domein
dx = B/(VB-1); dy = H/(VH-1); %Cell grotes
Cmet = 65; %Metaal
Cpla = 0.2; %Plastiek

%Stensil voor de inwendige punten 
%Zowel temperatuur en materiaal state zijn cell centered
%Verticale Faces

p = 5;
% MatArray = (1 - v) .^ p * Cpla + v .^ p * Cmet; %verschillende paper
MatDerArray = -p * (1 - v) .^(p-1) * Cpla + p * v .^ (p-1) * Cmet; %Andere paper

dg   = zeros(VB*VH, 1); % Diagonal
subdg  = zeros(VB*VH, 1); % Subdiagonal
supdg  = zeros(VB*VH, 1); % Superdiagonal
subvbdg = zeros(VB*VH, 1); % VB subdiagonal
supvbdg = zeros(VB*VH, 1); % VB superdiagoal

G = sparse(VB*VH,VB*VH);

for i = 1:VB-1 %in breedte
   %Onder
   k = i;
   
   dM1 = 1/2; dM2 = 1/2;
   R1 = MatDerArray(i,1)* dM1 * dy/dx/2;
   R2 = MatDerArray(i+1,1)* dM2 * dy/dx/2;
   
   G(k,   k) = G(k,   k) + R1 * SOL(k) - R1 * SOL(k+1); %Denk dat dit zou moeten
   G(k+1, k) = G(k+1, k) - R1 * SOL(k) + R1 * SOL(k+1);
    
   G(k,   k+1) = G(k,   k+1) + R2 * SOL(k) - R2 * SOL(k+1);
   G(k+1, k+1) = G(k+1, k+1) - R2 * SOL(k) + R2 * SOL(k+1);
   
   dg(k)   = dg(k) + R1 * SOL(k) - R1 * SOL(k+1);
   dg(k+1) = dg(k+1) - R2 * SOL(k) + R2 * SOL(k+1);
   subdg(k) =  subdg(k) - R1 * SOL(k) + R1 * SOL(k+1);
   supdg(k+1) = supdg(k+1)  + R2 * SOL(k) - R2 * SOL(k+1);
   
   %Intern
   for j = 2:VH-1 %in hooghte
        k = i+VB*(j-1);
        
        dM1 = 1/2; dM2 = 1/2;
        R1 = MatDerArray(i,j)* dM1 * dy/dx;
        R2 = MatDerArray(i+1,j)*dM2 * dy/dx;
        
        G(k,   k) = G(k,   k) + R1 * SOL(k) - R1 * SOL(k+1); %Denk dat dit zou moeten
        G(k+1, k) = G(k+1, k) - R1 * SOL(k) + R1 * SOL(k+1);
    
        G(k,   k+1) = G(k,   k+1) + R2 * SOL(k) - R2 * SOL(k+1);
        G(k+1, k+1) = G(k+1, k+1) - R2 * SOL(k) + R2 * SOL(k+1);
        
        dg(k)   = dg(k) + R1 * SOL(k) - R1 * SOL(k+1);
        dg(k+1) = dg(k+1) - R2 * SOL(k) + R2 * SOL(k+1);
        subdg(k) =  subdg(k) - R1 * SOL(k) + R1 * SOL(k+1);
        supdg(k+1) = supdg(k+1)  + R2 * SOL(k) - R2 * SOL(k+1);
   end
   
   %Boven
   k = i+VB*(VH-1);
   
   dM1 = 1/2; dM2 = 1/2;
   R1 = MatDerArray(i,VH)* dM1 * dy/dx/2;
   R2 = MatDerArray(i+1,VH)*dM2* dy/dx/2;
   
   G(k,   k) = G(k,   k) + R1 * SOL(k) - R1 * SOL(k+1); %Denk dat dit zou moeten
   G(k+1, k) = G(k+1, k) - R1 * SOL(k) + R1 * SOL(k+1);
    
   G(k,   k+1) = G(k,   k+1) + R2 * SOL(k) - R2 * SOL(k+1);
   G(k+1, k+1) = G(k+1, k+1) - R2 * SOL(k) + R2 * SOL(k+1);
   
   dg(k)   = dg(k) + R1 * SOL(k) - R1 * SOL(k+1);
   dg(k+1) = dg(k+1) - R2 * SOL(k) + R2 * SOL(k+1);
   subdg(k) =  subdg(k) - R1 * SOL(k) + R1 * SOL(k+1);
   supdg(k+1) = supdg(k+1)  + R2 * SOL(k) - R2 * SOL(k+1);
end

%Horizontale Faces
for j = 1:VH-1 %in hooghte
    %Links
    k = 1+VB*(j-1);
    
    dM1 = 1/2; dM2 = 1/2;
    R1 = MatDerArray(1,j)* dM1 * dy/dx/2;
    R2 = MatDerArray(1,j+1)*dM2 * dy/dx/2;
    
    G(k,   k) = G(k,   k) + R1 * SOL(k) - R1 * SOL(k+VB); %Denk dat dit zou moeten
    G(k+VB, k) = G(k+VB, k) - R1 * SOL(k) + R1 * SOL(k+VB);
    
    G(k,   k+VB) = G(k,   k+VB) + R2 * SOL(k) - R2 * SOL(k+VB);
    G(k+VB, k+VB) = G(k+VB, k+VB) - R2 * SOL(k) + R2 * SOL(k+VB);
    
    dg(k)   = dg(k) + R1 * SOL(k) - R1 * SOL(k+VB);
    dg(k+VB) = dg(k+VB) - R2 * SOL(k) + R2 * SOL(k+VB);
    subvbdg(k) =  subvbdg(k) - R1 * SOL(k) + R1 * SOL(k+VB);
    supvbdg(k+VB) = supvbdg(k+VB)  + R2 * SOL(k) - R2 * SOL(k+VB);
    
    %Intern
    for i = 2:VB-1 %in breedte
        k = i+VB*(j-1);

        dM1 = 1/2; dM2 = 1/2;
        R1 =  MatDerArray(i,j) * dM1 * dx/dy; 
        R2 = MatDerArray(i,j+1) * dM2 * dx/dy;
        
        G(k,   k) = G(k,   k) + R1 * SOL(k) - R1 * SOL(k+VB); 
        G(k+VB, k) = G(k+VB, k) - R1 * SOL(k) + R1 * SOL(k+VB);
    
        G(k,   k+VB) = G(k,   k+VB) + R2 * SOL(k) - R2 * SOL(k+VB);
        G(k+VB, k+VB) = G(k+VB, k+VB) - R2 * SOL(k) + R2 * SOL(k+VB);
        
        dg(k)   = dg(k) + R1 * SOL(k) - R1 * SOL(k+VB);
        dg(k+VB) = dg(k+VB) - R2 * SOL(k) + R2 * SOL(k+VB);
        subvbdg(k) =  subvbdg(k) - R1 * SOL(k) + R1 * SOL(k+VB);
        supvbdg(k+VB) = supvbdg(k+VB)  + R2 * SOL(k) - R2 * SOL(k+VB);
    end
   
    %Rechts
    k = VB+VB*(j-1);
    
    dM1 = 1/2; dM2 = 1/2;
    R1 = MatDerArray(VB,j)*dM1 * dx/dy/2; 
    R2 = MatDerArray(VB,j+1)* dM2 * dx/dy/2;
    
    G(k,   k) = G(k,   k) + R1 * SOL(k) - R1 * SOL(k+VB); 
    G(k+VB, k) = G(k+VB, k) - R1 * SOL(k) + R1 * SOL(k+VB);
    
    G(k,   k+VB) = G(k,   k+VB) + R2 * SOL(k) - R2 * SOL(k+VB);
    G(k+VB, k+VB) = G(k+VB, k+VB) - R2 * SOL(k) + R2 * SOL(k+VB);
    
    dg(k)   = dg(k) + R1 * SOL(k) - R1 * SOL(k+VB);
    dg(k+VB) = dg(k+VB) - R2 * SOL(k) + R2 * SOL(k+VB);
    subvbdg(k) =  subvbdg(k) - R1 * SOL(k) + R1 * SOL(k+VB);
    supvbdg(k+VB) = supvbdg(k+VB)  + R2 * SOL(k) - R2 * SOL(k+VB);
    
end
G = spdiags([subvbdg, subdg, dg, supdg, supvbdg], [-VB, -1, 0, 1, VB], VB*VH, VB*VH);
AG = L'*G;
end
