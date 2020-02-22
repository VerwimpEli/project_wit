function [f0val, df0dv, fval, dfdv] = heateq(v, M, VW, VH)
    % Solve heat equation and calculate gradient
    v = reshape(v, [VW, VH]);
    [Sol,K] = FVM(VW,VH, v);
    L = (K')\-ones(VW*VH,1);
    
    f0val = sum(Sol, 'all');   % TODO: Q!
    df0dv = Adjoint_Gradient(VW,VH,v,L,Sol)';

    % 1/n sum(v) <= M
    fval = sum(v, 'all') - M * VW * VH;  
    dfdv = ones(1, VW * VH);
end