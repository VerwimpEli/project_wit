function [f0val, df0dv, fval, dfdv] =  ...

    heateq(v, M, VW, VH, Q, Cmet, Cpla, BC0, BC1, BC2, BC3, p)

    if nargin < 12
        p = 1;
    end
                                                
    % Solve heat equation and calculate gradient
    v = reshape(v, [VW, VH]);
<<<<<<< HEAD:matlab/heateq.m
    [Sol,K] = FVM(VW,VH, v, Q, Cmet, Cpla, BC0, BC1, BC2, BC3, p);
    L = (K')\-scale(ones(VW*VH,1));
    
    f0val = sum(scale(Sol), 'all') ;   
    df0dv = Adjoint_Gradient(VW,VH,v,L,Sol, p, Cmet, Cpla)';
=======
    [Sol,K] = Harmonic_FVM(B,H,VW,VH, v, Q, Cmet, Cpla, BC0, BC1, BC2, BC3, p);
    L = (K')\-scale(ones(VW*VH,1));
    
    f0val = sum(scale(Sol), 'all') ;   
    df0dv = Harmonic_Adjoint_Gradient(B,H,VW,VH,v,L,Sol, p, Cmet, Cpla)';
>>>>>>> e1b68c7a123be878758721895d862934872543a0:matlab/Harmonic_heateq.m

    % 1/n sum(v) <= M
    fval = sum(scale(reshape(v, VW * VH, 1)), 'all') - M * VW * VH;  
    dfdv = scale(ones(VW * VH, 1))';
end