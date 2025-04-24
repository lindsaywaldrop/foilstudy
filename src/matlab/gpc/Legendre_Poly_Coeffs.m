
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: get Legendre Polynomial Coefficients: 
%
%       L_n(x) = c0 + c1*x + c2*x^2 +...+ c_n*x^n
%
%       input: n <- degree of Legendre poly
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function coeffs = Legendre_Poly_Coeffs(n)

    if n == 0

        coeffs = 1;      % L_0(x) = 1

    elseif n == 1

        coeffs = [1 0];  % L_1(x) = x

    else

        % Create Legendre Polynomials based off of recurrence relation  
        P_nm1 = 1;
        P_n = [1 0];
        for i=1:(n-1)
            P_np1 = ((2*i+1)*[P_n,0] - i*[0,0,P_nm1])/(i+1); % recurrence
            [P_nm1,P_n] = deal(P_n,P_np1); % shift
        end

        coeffs = P_np1;
        
    end