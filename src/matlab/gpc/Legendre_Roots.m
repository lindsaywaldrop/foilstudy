
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: get roots of (n+1)^(st) Legendre polynomial
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function leg_roots = Legendre_Roots(degree)

   % get coefficients of the nth degree Legendre polynomial
   coeffs = Legendre_Poly_Coeffs(degree);
   
   % find roots of polynomial with those coefficients
   leg_roots = roots(coeffs);