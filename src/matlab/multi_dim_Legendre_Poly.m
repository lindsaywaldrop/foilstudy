
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: compute multidimension Legendre Polynomial and evaluate for 
%               parameter combination in the vector paramVec
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = multi_dim_Legendre_Poly(alphaVec,paramVec)

    val = 1;
    for i=1:length(alphaVec)

        n = alphaVec(i);  % Get degree with specific single dimensional Leg. Poly  
        x = paramVec(i);  % Get value to evaluate single dimensional Leg. Poly at

        % Compute product of all Legendre Poly's evaluated at specific x and of
        %         particular degree
        val = val * Legendre_Poly(n,x);

    end