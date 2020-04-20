function AG = Adjoint_Gradient(VB,VH,v,L,SOL, p, Cmet, Cpla)%moet synchroon blijven met Bovenstaande Methode
H = 1; B = 1; %Hoogte en breedte van het domein
dx = B/(VB-1); dy = H/(VH-1); %Cell grotes

if nargin < 6  % Default p, Cmet and Cpla (backwards compatibility)
    p = 1;
    Cmet = 65;
    Cpla = 0.2;
end

%Stensil voor de inwendige punten 
%Zowel temperatuur en materiaal state zijn cell centered
%Verticale Faces

MatArray = (1 - v) .^ p * Cpla + v .^ p * Cmet; %verschillende paper
MatDerArray = -p * (1 - v) .^(p-1) * Cpla + p * v .^ (p-1) * Cmet; %Andere paper

dg   = zeros(VB*VH, 1); % Diagonal
subdg  = zeros(VB*VH, 1); % Subdiagonal
supdg  = zeros(VB*VH, 1); % Superdiagonal
subvbdg = zeros(VB*VH, 1); % VB subdiagonal
supvbdg = zeros(VB*VH, 1); % VB superdiagoal

for i = 1:VB-1 %in breedte
   %Onder
   k = i;
   M1 = MatArray(i,1); M2 = MatArray(i+1,1); %Materialen van de 2 cellen
   dM1 = 2/(M1^2*(1/M1+1/M2)^2); dM2 = 2/(M2^2*(1/M1+1/M2)^2);
   %(MatArray(i,1)+MatArray(i+1,1))/2;
   R1 = MatDerArray(i,1)* dM1 * dy/dx/2;
   R2 = MatDerArray(i+1,1)* dM2 * dy/dx/2;
   
   dg(k)   = dg(k) + R1 * SOL(k) - R1 * SOL(k+1);
   dg(k+1) = dg(k+1) - R2 * SOL(k) + R2 * SOL(k+1);
   subdg(k) =  subdg(k) - R1 * SOL(k) + R1 * SOL(k+1);
   supdg(k+1) = supdg(k+1)  + R2 * SOL(k) - R2 * SOL(k+1);
      
   %Intern
   for j = 2:VH-1 %in hooghte
        k = i+VB*(j-1);
        
        %M = (MatArray(i,j)+MatArray(i+1,j))/2;
        M1 = MatArray(i,j); M2 = MatArray(i+1,j); %Materialen van de 2 cellen
        dM1 = 2/(M1^2*(1/M1+1/M2)^2); dM2 = 2/(M2^2*(1/M1+1/M2)^2);
        R1 = MatDerArray(i,j)* dM1 * dy/dx;
        R2 = MatDerArray(i+1,j)*dM2 * dy/dx;
        
        dg(k)   = dg(k) + R1 * SOL(k) - R1 * SOL(k+1);
        dg(k+1) = dg(k+1) - R2 * SOL(k) + R2 * SOL(k+1);
        subdg(k) =  subdg(k) - R1 * SOL(k) + R1 * SOL(k+1);
        supdg(k+1) = supdg(k+1)  + R2 * SOL(k) - R2 * SOL(k+1);
   end
   
   %Boven
   k = i+VB*(VH-1);
   %M = (MatArray(i,VH)+MatArray(i+1,VH))/2;
   M1 = MatArray(i,VH); M2 = MatArray(i+1,VH); %Materialen van de 2 cellen
   dM1 = 2/(M1^2*(1/M1+1/M2)^2); dM2 = 2/(M2^2*(1/M1+1/M2)^2);
   R1 = MatDerArray(i,VH)* dM1 * dy/dx/2;
   R2 = MatDerArray(i+1,VH)*dM2* dy/dx/2;
   
    dg(k)   = dg(k) + R1 * SOL(k) - R1 * SOL(k+1);
    dg(k+1) = dg(k+1) - R2 * SOL(k) + R2 * SOL(k+1);
    subdg(k) =  subdg(k) - R1 * SOL(k) + R1 * SOL(k+1);
    supdg(k+1) = supdg(k+1)  + R2 * SOL(k) - R2 * SOL(k+1);
   
end

%Horizontale Faces
for j = 1:VH-1 %in hooghte
    %Links
    k = 1+VB*(j-1);
    
    %M = (MatArray(1,j)+MatArray(1,j+1))/2; van algebraisch gemiddelde
    M1 = MatArray(1,j); M2 = MatArray(1,j+1); %Materialen van de 2 cellen
    dM1 = 2/(M1^2*(1/M1+1/M2)^2); dM2 = 2/(M2^2*(1/M1+1/M2)^2);
    R1 = MatDerArray(1,j)* dM1 * dy/dx/2;
    R2 = MatDerArray(1,j+1)*dM2 * dy/dx/2;
    
    dg(k)   = dg(k) + R1 * SOL(k) - R1 * SOL(k+VB);
    dg(k+VB) = dg(k+VB) - R2 * SOL(k) + R2 * SOL(k+VB);
    subvbdg(k) =  subvbdg(k) - R1 * SOL(k) + R1 * SOL(k+VB);
    supvbdg(k+VB) = supvbdg(k+VB)  + R2 * SOL(k) - R2 * SOL(k+VB);
       
    %Intern
    for i = 2:VB-1 %in breedte
        k = i+VB*(j-1);
        %M = (MatArray(i,j)+MatArray(i,j+1))/2;
        M1 = MatArray(i,j); M2 = MatArray(i,j+1);
        dM1 = 2/(M1^2*(1/M1+1/M2)^2); dM2 = 2/(M2^2*(1/M1+1/M2)^2);
        R1 =  MatDerArray(i,j) * dM1 * dx/dy; 
        R2 = MatDerArray(i,j+1) * dM2 * dx/dy;
        
        dg(k)   = dg(k) + R1 * SOL(k) - R1 * SOL(k+VB);
        dg(k+VB) = dg(k+VB) - R2 * SOL(k) + R2 * SOL(k+VB);
        subvbdg(k) =  subvbdg(k) - R1 * SOL(k) + R1 * SOL(k+VB);
        supvbdg(k+VB) = supvbdg(k+VB)  + R2 * SOL(k) - R2 * SOL(k+VB);
    end
   
    %Rechts
    k = VB+VB*(j-1);
    
    %M = (MatArray(VB,j)+MatArray(VB,j+1))/2;
    M1 = MatArray(VB,j); M2 = MatArray(VB,j+1);
    dM1 = 2/(M1^2*(1/M1+1/M2)^2); dM2 = 2/(M2^2*(1/M1+1/M2)^2);
    R1 = MatDerArray(VB,j)*dM1 * dx/dy/2; 
    R2 = MatDerArray(VB,j+1)* dM2 * dx/dy/2;
    
    dg(k)   = dg(k) + R1 * SOL(k) - R1 * SOL(k+VB);
    dg(k+VB) = dg(k+VB) - R2 * SOL(k) + R2 * SOL(k+VB);
    subvbdg(k) =  subvbdg(k) - R1 * SOL(k) + R1 * SOL(k+VB);
    supvbdg(k+VB) = supvbdg(k+VB)  + R2 * SOL(k) - R2 * SOL(k+VB);
    
end
G = spdiags([subvbdg, subdg, dg, supdg, supvbdg], [-VB, -1, 0, 1, VB], VB*VH, VB*VH);
AG = L'*G;
end
