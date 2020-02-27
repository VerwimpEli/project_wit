% Initalize problem
VW = 50;
VH = 50; 
Q  = 200;
Cpla = 0.2;
Cmet = 65;
M  = 0.4;   % Metal to plastic ratio
v  = rand(VW * VH, 1) * M;


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
maxiter = 200;
kkttol  = 1e-8;
kktnorm = 1.0;

% Inital 
v0 = v;
[f0val,f0grad,fval,fgrad] = heateq(v, M, VW, VH, Q, Cmet, Cpla);
fvals = zeros(maxiter, 1);
symm  = zeros(maxiter, 1);
fvals(1) = f0val;
symm(1)     = norm(reshape(v0, VW, VH) - flip(reshape(v0, VW, VH)), 'fro');


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
    
    [f0val,f0grad,fval,fgrad] = heateq(v, M, VW, VH, Q, Cmet, Cpla);
    fvals(iter) = f0val;
    symm(iter)     = norm(reshape(v, VW, VH) - flip(reshape(v, VW, VH)), 'fro');
    
    
    % Check KKT conditions for optimality
    [residu,kktnorm,residumax] = ...
    kktcheck(m, n, vmma, ymma, zmma, lam, xsi, eta, mu, zet, s, ...
             vmin, vmax, f0grad, fval, fgrad, a0, a, c, d);
end



% --- Plots ----


v = reshape(v, VW, VH);
figure
imagesc(v)
colorbar;

% Inital temp
tsol0 = FVM(VW, VH, reshape(v0, VW, VH), Q, Cmet, Cpla);
% Final temp
tsol  = FVM(VW, VH, v, Q, Cmet, Cpla);

maxt = max(tsol0);
mint = min(tsol0);

figure
subplot(121);
imagesc(reshape(tsol0, VW, VH));
caxis([mint, maxt]);
grid on

subplot(122);
imagesc(reshape(tsol, VW, VH));
caxis([mint, maxt]);
grid on

colorbar;







