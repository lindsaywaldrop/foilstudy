%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: DERIVATIVE of the specified activation function
%           
%       NOTE: [1] setup as function of matrix (or vector)
%             [2] make sure activation_function.m uses same
%                 activation function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = act_function_PRIME( Z )

%--------------------------------
% SIGMOID DERIVATIVE
%--------------------------------
%A = exp( -Z ) ./ ( 1 + exp( -Z ) ).^2;

%--------------------------------
% ReLU DERIVATIVE
%--------------------------------
A = zeros(size(Z));
%
A( find(Z>0) ) = 1;