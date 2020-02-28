function AG = Adjoint_Gradient(VB,VH,v,L,SOL)%moet synchroon blijven met Bovenstaande Methode
H = 1; B = 1; %Hoogte en breedte van het domein
dx = B/(VB-1); dy = H/(VH-1); %Cell grotes
Cmet = 65; %Metaal
Cpla = 0.2; %Plastiek

% G = zeros(VB*VH,VB*VH);
dg   = zeros(VB*VH, 1); % Diagonal
sdg  = zeros(VB*VH, 1); % Subdiagonal
vbdg = zeros(VB*VH, 1); % VB diagonal

%Stensil voor de inwendige punten 
%Zowel temperatuur en materiaal state zijn cell centered
%Verticale Faces

p = 5;
MatArray = p * (1 - v) .^(p-1) * Cpla + p * v .^ (p-1) * Cmet;

for i = 1:VB-1 %in breedte
   %Onder
   k = i;
   M = (MatArray(i,1)+MatArray(i+1,1))/2;
   R1 = M * dy/dx/4;
   R2 = M * dy/dx/4;
   
%    G(k,   k) = G(k,   k) + R1 * SOL(k) - R2 * SOL(k+1);
%    G(k+1, k) = G(k+1, k) - R1 * SOL(k) + R2 * SOL(k+1);
%    
%    G(k,   k+1) = G(k,   k+1) + R1 * SOL(k) - R2 * SOL(k+1);
%    G(k+1, k+1) = G(k+1, k+1) - R1 * SOL(k) + R2 * SOL(k+1);
   
   dg(k)   = dg(k) + R1 * SOL(k) - R2 * SOL(k+1);
   dg(k+1) = dg(k+1) - R1 * SOL(k) + R2 * SOL(k+1);
   sdg(k)  = sdg(k) -  R1 * SOL(k) + R2 * SOL(k+1);
   
   %Intern
   for j = 2:VH-1 %in hooghte
        k = i+VB*(j-1);
        
        M = (MatArray(i,j)+MatArray(i+1,j))/2;
        R1 = M * dy/dx/4*(Cmet-Cpla);
        R2 = M * dy/dx/4*(Cmet-Cpla);
       
       dg(k)   = dg(k) + R1 * SOL(k) - R2 * SOL(k+1);
       dg(k+1) = dg(k+1) - R1 * SOL(k) + R2 * SOL(k+1);
       sdg(k)  = sdg(k) -  R1 * SOL(k) + R2 * SOL(k+1);
   end
   
   %Boven
   k = i+VB*(VH-1);
   
   M = (MatArray(i,VH)+MatArray(i+1,VH))/2;
   R1 = M * dy/dx/4*(Cmet-Cpla);
   R2 = M * dy/dx/4*(Cmet-Cpla);
   
   dg(k)   = dg(k) + R1 * SOL(k) - R2 * SOL(k+1);
   dg(k+1) = dg(k+1) - R1 * SOL(k) + R2 * SOL(k+1);
   sdg(k)  = sdg(k) -  R1 * SOL(k) + R2 * SOL(k+1);
end

%Horizontale Faces
for j = 1:VH-1 %in hooghte
    %Links
    k = 1+VB*(j-1);
    
    M = (MatArray(1,j)+MatArray(1,j+1))/2;
    R1 = M * dy/dx/4*(Cmet-Cpla);
    R2 = M * dy/dx/4*(Cmet-Cpla);
    
    dg(k)    = dg(k) + R1 * SOL(k) - R2 * SOL(k+VB);
    dg(k+VB) = dg(k+VB) - R1 * SOL(k) + R2 * SOL(k+VB);
    vbdg(k)  = vbdg(k) -  R1 * SOL(k) + R2 * SOL(k+VB);
    
    %Intern
    for i = 2:VB-1 %in breedte
        k = i+VB*(j-1);
        M = (MatArray(i,j)+MatArray(i,j+1))/2;
        R1 = M * dx/dy/2*(Cmet-Cpla); 
        R2 = M * dx/dy/2*(Cmet-Cpla);
        
        dg(k)    = dg(k) + R1 * SOL(k) - R2 * SOL(k+VB);
        dg(k+VB) = dg(k+VB) - R1 * SOL(k) + R2 * SOL(k+VB);
        vbdg(k)  = vbdg(k) -  R1 * SOL(k) + R2 * SOL(k+VB);
    end
   
    %Rechts
    k = VB+VB*(j-1);
    
    M = (MatArray(VB,j)+MatArray(VB,j+1))/2;
    R1 = M * dx/dy/4*(Cmet-Cpla); 
    R2 = M * dx/dy/4*(Cmet-Cpla);
    
    dg(k)    = dg(k) + R1 * SOL(k) - R2 * SOL(k+VB);
    dg(k+VB) = dg(k+VB) - R1 * SOL(k) + R2 * SOL(k+VB);
    vbdg(k)  = vbdg(k) -  R1 * SOL(k) + R2 * SOL(k+VB);
end

G = spdiags([vbdg, sdg, dg, [0; -sdg(1:end-1)], ...
            [zeros(VB, 1); -vbdg(1:end-VB)]], ...
            [-VB, -1, 0, 1, VB], VB*VH, VB*VH);
AG = L'*G;
end
