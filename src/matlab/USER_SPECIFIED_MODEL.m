%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computes: (1) Performs a model evaulation for a particular parameter
%               combination, in the input VEC
%
%           (2) Since Legendre Polys are defined in [-1,1], transform the 
%               VEC input from [-1,1] to [A,B] accordingly for model
%               evaluation
%
%     Author: Nick Battista
%     Date: 11/11/2021
%     Institution: TCNJ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function val = USER_SPECIFIED_MODEL( VEC )

%-------------------------------------------------------------------------
% MODEL: polynomial model (ex1 from Sudret 2008)
%
%       Y = 1/2^N * PROD_i ( 3*x_i^2 + 1 )
%
%       NOTE: 1. x_i in [-1,1]
%             2. Validated against Sudret for 3 parameters
%
%       Sobol Indices (for 3 parameters):
%          -> 1st Order: all 0.274725274725275
%          -> 2nd Order: all 0.0549450549450549
%          -> 3rd Order: 0.010989010989011
%
%-------------------------------------------------------------------------
% N = length(VEC);
% prod = 1;
% for i=1:N
%     prod = prod * ( 3*( VEC(i) )^2 + 1 );
% end
% val =  1/2^N * prod;


%-------------------------------------------------------------------------
% MODEL: Ishigami function (ex2 from Sudret 2008)
%
%       Y = sin(x1) + a*sin^2(x2) + b*x3^4*sin(x1)
%
%       NOTE: 1. a,b model parameters (not varied)
%             2. x1,x2,x3 in [-pi,pi]
%             2. Validated against Sudret for 3 parameters
%
%       ACTUAL SOBOL INDICES (for 3 parameters):
%          -> 1ST-ORDER INDICES:  [0.3138, 0.4424, 0]
%          -> 2ND-ORDER INDICES:  [0 0.2436 0]
%          -> s123-INDEX:         [0]
%          -> TOTAL-ORDER INDICES:[0.5574, 0.4424, 0.2436] 
%
%-------------------------------------------------------------------------

% Define set parameters of Ishigami Function
a = 7;
b = 0.1;

% Transform VEC values from [-1,1] -> [-pi,pi] = [A,B]
A = -pi;
B = pi;  
VEC_T = transform_Values_To_Desired_Interval(VEC,A,B); % Get Transformed Evaulation Pts

% Evaluate Ishigami function
val = sin(VEC_T(1)) + a*( sin(VEC_T(2)) )^2 + b*(VEC_T(3))^4*sin(VEC_T(1));



%-------------------------------------------------------------------------
% MODEL: Sobol function (ex3 from Sudret 2008)
%
%       Y = PROD_i=1^N ( abs(4*x_i - 2) + a_i ) / (1 + a_i)
%
%       a = [1 2 5 10 20 50 100 500]
%       x_i \in [0,1]
%
%       ACTUAL SOBOL INDICES (for 3 parameters):
%          -> 1ST-ORDER INDICES:  [0.6037 0.2683 0.0671 0.0200 0.0055 0.0009 0.0002 0.0000]
%          -> TOTAL-ORDER INDICES:[0.6342 0.2945 0.0756 0.0227 0.0062 0.0011 0.0003 0.0000]
%             (see Sudret 2008 for other indices)
%
%      **** FOR THIS EXAMPLE, USE MORE TRAINING DATA SO  ****
%      **** INFORMATION MATRIX ISN'T SINGULAR!!!         ****
%      ****     ex) (p,N_subset)=(3, 3*(N-1)*cap_P ),    ****
%      ****         (p,N_subset)=(4, 3*(N-1)*cap_P )     ****
%             
%-------------------------------------------------------------------------
% 
% % Sobol' function parameters (from Sudret 2008)
% aVec = [1 2 5 10 20 50 100 500];
% 
% % Transform VEC values from [-1,1] -> [0,1] = [A,B]
% A = 0;
% B = 1;  
% VEC_T = transform_Values_To_Desired_Interval(VEC,A,B); % Get Transformed Evaulation Pts
% 
% % Sobol' function evaluation
% prod=1;
% for n=1:length(VEC_T)
%     prod = prod * ( abs( 4*VEC_T(n) - 2) + aVec(n) ) / (1 + aVec(n) );
% end
% val = prod;




%-------------------------------------------------------------------------
% MODEL: HIGHLY Oscillatory Ishigami function 
%
%       M = 2,3,4... to test degree, d, needed to capture accurate
%                   response surface
%
%       Y = sin(x1) + a*sin^2( M*x2 ) + b*x3^4*sin(x1)
%
%       NOTE: 1. a,b model parameters (not varied)
%             2. x1,x2,x3 in [-pi,pi]
%             3. Example meant for 3 parameters
%             4. Example to show when higher-order polynomials
%                are necessary for the surrogate 
%
%-------------------------------------------------------------------------
% 
% % Define set parameters of Ishigami Function
% a = 7;
% b = 0.2;
% M = 2.6;
% 
% % Transform VEC values from [-1,1] -> [-pi,pi] = [A,B]
% A = -pi;
% B = pi;  
% VEC_T = transform_Values_To_Desired_Interval(VEC,A,B); % Get Transformed Evaulation Pts
% 
% %Evaluate HIGHLY OSCILLATORY Ishigami function
% val = VEC_T(1) + a*( sin(M*VEC_T(2)-0.351) )^2 + b*(VEC_T(3))^2*VEC_T(1)^2 + 5;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: transforms values in [-1,1] to desired interval [A,B]
%           (when model function inputs defined over a different
%                                               domain than [-1,1] )
%
%       TRANSFORM: y=mx+b
%           y(-1) = -m + b = A
%           y(1) = m + b = B
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function VEC_T = transform_Values_To_Desired_Interval(VEC,A,B)


    % Setup linear system
    MAT = [-1 1; 1 1];
    RHS = [A;B];
    coeffs = MAT \ RHS;

    m = coeffs(1);
    b = coeffs(2);

    VEC_T = m*VEC + b;



