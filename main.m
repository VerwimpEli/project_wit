% Initalize problem
VW = 25;
VH = 25; 
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
maxiter = 10;
kkttol  = 1e-8;
kktnorm = 1.0;

% Inital 
v0 = v;
[f0val,f0grad,fval,fgrad] = heateq(v, M, VW, VH);

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
    
    [f0val,f0grad,fval,fgrad] = heateq(v, M, VW, VH);
    
    % Check KKT conditions for optimality
    [residu,kktnorm,residumax] = ...
    kktcheck(m, n, vmma, ymma, zmma, lam, xsi, eta, mu, zet, s, ...
             vmin, vmax, f0grad, fval, fgrad, a0, a, c, d);
end



% --- Plots ----


v = reshape(v, VW, VH);
figure
surf(v)

% Inital temp
tsol0 = FVM(VW, VH, reshape(v0, VW, VH));
% Final temp
tsol  = FVM(VW, VH, v);

figure
subplot(121);
surf(reshape(tsol0, VW, VH));
caxis([0, 2.05]);

subplot(122);
surf(reshape(tsol, VW, VH));
caxis([0, 2.05]);

colorbar;







