
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Evaluate Multi-dimensional Legendre Polynomial at specific
%           x,y,z
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = evaluate_MultiDim_Legendre_Poly(sCoeffs,alphaMAT,VEC)

    sum = 0;
    for i=1:length(alphaMAT(:,1))

        prod = 1;
        for j=1:length(alphaMAT(1,:))
            prod = prod * Legendre_Poly(alphaMAT(i,j),VEC(j)); 
        end

        sum = sum + sCoeffs(i)*prod;

    end

    val = sum;
