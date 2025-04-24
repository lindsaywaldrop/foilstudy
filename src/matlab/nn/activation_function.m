%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Specified Activation function
%           
%       NOTE: [1] setup as function of matrix (or vector)
%             [2] make sure act_function_PRIME.m uses same 
%                 activation function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = activation_function( Z )

%--------------------------------
% SIGMOID ACTIVATION
%--------------------------------
%A = 1 ./ ( 1 + exp( -Z ) );

%--------------------------------
% ReLU
%--------------------------------
A=max(0,Z);

