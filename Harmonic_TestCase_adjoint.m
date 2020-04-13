clear;
%Testcase for the testing the adjoint method (used for calculation of the
%gradient) by comparision with a finite difference approxition. Also used
%to generate a FVM solution and gradient for comparison with the cpp
%implementation
%REMARK 1 check if the cpp programme is using the same parameters
%Q,BC's,...
%REMARK 2 Use a non -random, non homogeneous material, to check if all
%matrial selection indices are correct.
%REMARK 3 Location which cells belong to a given BC is implemented
%different is cpp. Therefore some difference might be possible for special meshes. 
%This is not the case for (eg 25x25)
VH = 25; VW = 25; % Aantal volumes in de hoogte en breedte. Incluisief de kleinere op de randen
Varray = linspace(0.4,0.6,VW*VH)'; %Equispaced tussen 0.4 en 0.6
v = reshape(Varray, [VW, VH]);

%%%Boundary condition
VDB = 0.3 * VH + 1;   % Dirichlet begin and end
VDE = VH - 0.3 * VH;    

BC0 = [['N',1,1,0]; ['N',2,VW-1,0]; ['N',VW,VW,0]];     % Onder
BC1 = [['N',1,1,0]; ['N',2,VDB-1,0]; ['D',VDB, VDE, 293]; ['N',VDE+1,VH-1, 0]; ['N',VH,VH,0]];    % Rechts 
BC2 = [['N',1,1,0]; ['N',2,VW-1,0]; ['N',VW,VW,0]];     % Boven
BC3 = [['N',1,1,0]; ['N',2,VDB-1,0]; ['D',VDB, VDE, 293]; ['N',VDE+1,VH-1, 0]; ['N',VH,VH,0]];    % Links 
%%Extra parameters
q = 200;
Cmet = 65; Cpla = 0.2;
p = 2;

%Refentie oplossing
[Sol,K] = Harmonic_FVM(VW, VH, v, q, Cmet, Cpla, BC0, BC1, BC2, BC3, p);


%Visualisatie %Heel ruw komt niet direct overeen met echte systeem door
%speciale Mesh
SOL = reshape(Sol,[VW,VH]);
% figure(1);
% surf(SOL); 
% xlabel("X"); ylabel("Y"); zlabel("Temperatuur")

%Benaderen Jacobiaan
FDG = FD_G(VW,VH,Varray,q,Cmet, Cpla, BC0, BC1, BC2, BC3, p)';

%Adjoint
L = (K')\-scale(ones(VW*VH,1));
AG = Harmonic_Adjoint_Gradient(VW,VH,v,L,Sol, p, Cmet, Cpla);

ERR1 = log10(abs(reshape(AG-FDG,[VW,VH])./reshape(AG+FDG,[VW,VH])/2));


figure(2);
subplot(1,3,1); 
surf(reshape(AG,[VW,VH])); 
title("AG");
subplot(1,3,2)
surf(ERR1); 
title("10 LOG REL ERR1 TOV AVG");
xlabel("X"); ylabel("Y"); zlabel("ERROR1")
subplot(1,3,3)
surf(reshape(FDG,[VW,VH])); 
title("FDG");

%%%Bereken van de gradient DMV Finite difference
function J = FD_G(VB,VH,Varray,q,Cmet, Cpla, BC0, BC1, BC2, BC3, p)
    Delta = 10^-6;
    J = zeros(size(Varray));
    %Referentie Oplossing
    [Sol, ~, ~, ~] =  Harmonic_heateq(Varray, 1, VB, VH, q, Cmet, Cpla, BC0, BC1, BC2, BC3, p);
    
    %Loop
    for i = 1:size(Varray,1)
            Varray2 = Varray; Varray2(i) = Varray2(i)*(1+Delta);
            [FD_Sol, ~, ~, ~] =  Harmonic_heateq(Varray2, 1, VB, VH, q, Cmet, Cpla, BC0, BC1, BC2, BC3, p);
            J(i)= (FD_Sol - Sol)/(Delta*Varray(i)); %of is dit de gradient?
    end
end


%%%Gradient in materiaal Array
function v = material_gradient(VW,VH,Start,End)
    v = ones(VW,VH);
    for i = 1:VW
        v(:,i) = Start + (i-1)/(VW-1)*(End-Start);
    end
end