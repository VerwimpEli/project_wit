function [f0val, df0dv, fval, dfdv] =  ...
    Harmonic_heateq(v, M, VW, VH, Q, Cmet, Cpla, BC0, BC1, BC2, BC3, p)
                                                
    % Solve heat equation and calculate gradient
    v = reshape(v, [VW, VH]);
    [Sol,K] = Harmonic_FVM(VW,VH, v, Q, Cmet, Cpla, BC0, BC1, BC2, BC3, p);
    L = (K')\-scale(ones(VW*VH,1));
    
    f0val = sum(scale(Sol), 'all') ;   
    df0dv = Harmonic_Adjoint_Gradient(VW,VH,v,L,Sol, p, Cmet, Cpla)';

    % 1/n sum(v) <= M
    fval = sum(scale(reshape(v, VW * VH, 1)), 'all') - M * VW * VH;  
    dfdv = scale(ones(VW * VH, 1))';
end