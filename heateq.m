function [f0val, df0dv, fval, dfdv] = heateq(v, M, VW, VH, Q, Cmet, Cpla)
    % Solve heat equation and calculate gradient
    v = reshape(v, [VW, VH]);
    [Sol,K] = FVM(VW,VH, v, Q, Cmet, Cpla);
    L = (K')\-scale(ones(VW*VH,1));
    
    factor = 2;
    v = reshape(v, VW * VH, 1);
    f0val = sum(scale(Sol) - factor * v .* (v - 1), 'all') ;   % TODO: Q!
    df0dv = Adjoint_Gradient(VW,VH,v,L,Sol)' - 2 * factor * v - factor;

    % 1/n sum(v) <= M
    fval = sum(scale(reshape(v, VW * VH, 1)), 'all') - M * VW * VH;  
    dfdv = ones(1, VW * VH);
end