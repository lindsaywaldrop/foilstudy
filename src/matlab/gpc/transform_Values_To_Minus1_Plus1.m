function VAL_T = transform_Values_To_Minus1_Plus1(VAL,A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: transforms values in [A,B] to desired interval [-1,1]
%
%       TRANSFORM: y=mx+b
%           y(A) = m(A) + b = -1
%           y(B) = m(B) + b = 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Setup linear system
MAT = [A 1; B 1];
RHS = [-1;1];
coeffs = MAT \ RHS;

m = coeffs(1);
b = coeffs(2);

VAL_T = m*VAL + b;
