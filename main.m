% Initalize problem
VW = 128;
VH = 128; 
Q  = 2000;
Cpla = 0.2;
Cmet = 65;
M  = 0.2;   % Metal to plastic ratio
% v  = rand(VW * VH, 1) * M;
p = 1;
pmax = 5;
v = single(zeros(VW * VH, 1));

% Boundary Conditions
VDB = 0.3 * VH + 1;   % Dirichlet begin and end
VDE = VH - 0.3 * VH;

% Assignment 

BC0 = [['N',1,1,0]; ['N',2,VW-1,0]; ['N',VW,VW,0]];     % Onder
BC1 = [['N',1,1,0]; ['N',2,VDB-1,0]; ['D',VDB, VDE, 293]; ['N',VDE+1,VH-1, 0]; ['N',VH,VH,0]];    % Rechts 
BC2 = [['N',1,1,0]; ['N',2,VW-1,0]; ['N',VW,VW,0]];     % Boven
BC3 = [['N',1,1,0]; ['N',2,VDB-1,0]; ['D',VDB, VDE, 293]; ['N',VDE+1,VH-1, 0]; ['N',VH,VH,0]];    % Links 

% Finite Volume paper
% 
% BC0 = [['N',1,1,0]; ['N',2,VW-1,0]; ['N',VW,VW,0]];     % Onder
% BC1 = [['N',1,1,0]; ['N',2,VW-1,0]; ['N',VH,VH,0]];    % Rechts 
% BC2 = [['D',1,1,273]; ['D',2,VW-1,273]; ['D',VW,VW,273]];     % Boven
% BC3 = [['D',1,1,273]; ['D',2,VW-1,273]; ['D',VW,VW,273]];    % Links 

% BC0 = [['N',1,1,0]; ['N',2,VW-1,0]; ['N',VW,VW,0]];     % Onder
% BC1 = [['N',1,1,0]; ['N',2,VW-1,0]; ['N',VH,VH,0]];    % Rechts 
% BC2 = [['N',1,1,0]; ['N',2,VW-1,0]; ['N',VW,VW,0]];     % Boven
% sw = 2;
% BC3 = [['N',1,1,0]; ['N',2,VW/2-sw-1,0]; ['D',VW/2-sw,VW/2+sw,293]; ...
%        ['N',VW/2+sw+1, VW-1, 0]; ['N',VW,VW,0]];    % Links 

% Initialize MMA optimizer
m       = 1;
n       = VW * VH;
epsimin = 0.0000001;
vold1   = v;
vold2   = v;
vmin    = zeros(n, 1);
vmax    = ones(n, 1);
low     = vmin;
upp     = vmax;
c       = 1e8;
d       = 1;
a0      = 1;
a       = 0;
iter    = 0;
maxiter = 10;
kkttol  = 1e-8;
kktnorm = 1.0;

% Inital 
v0 = v;
[f0val,f0grad,fval,fgrad] = Harmonic_heateq(v, M, VW, VH, Q, Cmet, Cpla, BC0, BC1, BC2, BC3, p);
fvals = zeros(maxiter, 1);
symm  = zeros(maxiter, 1);
fvals(1) = f0val;
symm(1)  = norm(reshape(v0, VW, VH) - flip(reshape(v0, VW, VH)), 'fro');


% Loop
while kktnorm > kkttol && iter < maxiter 
    iter = iter + 1;
    disp("Iteration: " + iter);
    % mma optimizer
    [vmma, ymma, zmma, lam, xsi, eta, mu, zet, s, low, upp] = ...
        mmasub(m, n, iter, v, vmin, vmax, vold1, vold2, ...
               f0val, f0grad, fval, fgrad, low, upp, a0, a, c, d);
    
    vold2 = vold1;
    vold1 = v;
    v = vmma;
    
    [f0val,f0grad,fval,fgrad] = Harmonic_heateq(v, M, VW, VH, Q, Cmet, Cpla, BC0, BC1, BC2, BC3, p);
    fvals(iter+1) = f0val;
    symm(iter+1)  = norm(reshape(v, VW, VH) - flip(reshape(v, VW, VH)), 'fro');
    
    
    % Check KKT conditions for optimality
    [residu,kktnorm,residumax] = ...
    kktcheck(m, n, vmma, ymma, zmma, lam, xsi, eta, mu, zet, s, ...
             vmin, vmax, f0grad, fval, fgrad, a0, a, c, d);
         
    if abs(fvals(iter+1) - fvals(iter)) / (VW * VH) < 0.01 && p < pmax 
        p = p + 1;
        disp("P increased: " + p)
    end
end



% --- Plots ----


v = reshape(v, VW, VH);
figure
imagesc(flip(v'))
colorbar;

% Inital temp
tsol0 = FVM(VW, VH, reshape(v0, VW, VH), Q, Cmet, Cpla, BC0, BC1, BC2, BC3);
% Final temp
tsol  = FVM(VW, VH, v, Q, Cmet, Cpla, BC0, BC1, BC2, BC3);

maxt = max(tsol0);
mint = min(tsol0);

figure
subplot(121);
imagesc(flip(reshape(tsol0, VW, VH)'));
caxis([273, maxt]);
grid on

subplot(122);
imagesc(flip(reshape(tsol, VW, VH)'));
caxis([273, maxt]);
grid on
colorbar;

figure
subplot(121);
plot(fvals);
subplot(122);
plot(symm);








