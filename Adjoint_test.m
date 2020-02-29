% Initalize problem
VW = 10;
VH = 10; 
Q  = 2000;
Cpla = 0.2;
Cmet = 65;
M  = 0.2;   % Metal to plastic ratio
v  = rand(VW * VH, 1) * M;
% v = zeros(VW * VH, 1);

% Boundary Conditions
VDB = 0.3 * VH + 1;   % Dirichlet begin and end
VDE = VH - 0.3 * VH;

% Heat equation solved with Adjoint mehod

BC0 = [['N',1,1,0]; ['N',2,VW-1,0]; ['N',VW,VW,0]];     % Onder
BC1 = [['N',1,1,0]; ['N',2,VDB-1,0]; ['D',VDB, VDE, 293]; ['N',VDE+1,VH-1, 0]; ['N',VH,VH,0]];    % Rechts 
BC2 = [['N',1,1,0]; ['N',2,VW-1,0]; ['N',VW,VW,0]];     % Boven
BC3 = [['N',1,1,0]; ['N',2,VDB-1,0]; ['D',VDB, VDE, 293]; ['N',VDE+1,VH-1, 0]; ['N',VH,VH,0]];    % Links 

[f0val,f0grad,fval,fgrad] = heateq(v, M, VW, VH, Q, Cmet, Cpla, BC0, BC1, BC2, BC3);

% Heat equation solved with finite differences

J = FD_J(v, M, VW, VH, Q, Cmet, Cpla, BC0, BC1, BC2, BC3);
G_FDT = J;

err = norm(G_FDT-f0grad);

ERR1 = G_FDT-f0grad;
rel_err = reshape(ERR1./f0grad,[VW,VH]);
% rel_err = reshape(ERR1,[VW,VH]);

figure(1);
subplot(1,3,1);
surf(reshape(f0grad,[VW,VH])); 
xlabel("X"); ylabel("Y"); zlabel("gradient")
title("Adjoint");
subplot(1,3,2)
surf(reshape(G_FDT,[VW,VH])); 
title("Finite differences");
xlabel("X"); ylabel("Y"); zlabel("gradient")
subplot(1,3,3)

surf(rel_err);
title("Relative Error");
xlabel("X"); ylabel("Y"); zlabel("error");

figure(2)
surf(log10(abs(rel_err)));
% semilogy(sort(reshape(ERR1./f0grad,[VW,VH])));


% Calculate gradient of cost function using finite differences
function J = FD_J(v, M, VW, VH, Q, Cmet, Cpla, BC0, BC1, BC2, BC3)
    Delta = 10^-6;
    J = zeros(size(v));
    [f0val1,f0grad,fval,fgrad] = heateq(v, M, VW, VH, Q, Cmet, Cpla, BC0, BC1, BC2, BC3);
    %Loop
    for i = 1:size(v,1)
        Varray2 = v;
        Varray2(i) = Varray2(i)*(1+Delta);
        [f0val2,f0grad,fval,fgrad] = heateq(Varray2, M, VW, VH, Q, Cmet, Cpla, BC0, BC1, BC2, BC3);
        J(i)= (f0val2 - f0val1)/(Delta*v(i)); %of is dit de gradient?
    end
end

