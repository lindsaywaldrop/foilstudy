%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Computes: (1) the Generalized Polynomial Chaos (gPC) expansion of the
%               model function (from the USER_SPECIFIED_MODEL) script.
%
%           (2) Plots the response surface for a particular 2D subspace
%
%           (3) From the gPC expansion, the Sobol sensitivity indices are
%               calculated from the gPC coefficients
%
%     Author: Nick Battista
%     Date: 11/11/2021
%     Institution: TCNJ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [param_combo_full, param_combo] = USER_RUN_gPC_Expansion()


%-------------------------------------------------------
% SAVE gPC Data (to '.mat' file)
%-------------------------------------------------------
save_Info = 0;

%-------------------------------------------------------
% TEST gPC Expansion 
%  (on specific 2D subspace, see test_gPC_Expansion.m)
%-------------------------------------------------------
test_Flag = 1;


%-------------------------------------------------------
% SET gPC PARAMETERS
% (this code only setup for 3 parameters)
%-------------------------------------------------------
N = 3;                   % 3 uncertain parameters
p = 6;                   % highest degree multivariable polynomial power
cap_P = nchoosek(N+p,p); % total # of terms in gPC expansion


%-------------------------------------------------------
% FIND COLLOCATION POINTS
%       (roots of (p+1)st degree Legendre polynomial)
%-------------------------------------------------------
poly_roots = Legendre_Roots(p+1);


%-------------------------------------------------------
% FIND COMBINATIONS OF ALL COLLOCATIONS PTS
%-------------------------------------------------------
param_combo_full = compute_All_Collo_PT_Combinations(N,poly_roots);


%------------------------------------------------------------
% GET SUBSET OF ALL PARAMETER COMBINATIONS FOR LEAST SQUARES
%------------------------------------------------------------
N_subset = 1*(N-1)*cap_P;
N_full = length(param_combo_full);
%param_combo = sample_Parameter_Combinations(N_subset,param_combo_full);
param_combo = param_combo_full;


%------------------------------------------------------------
% PRINT SOME GENERAL INFORMATION TO SCREEN
%------------------------------------------------------------
fprintf('\n\n--------------------------------------\n');
fprintf('gPC Expansion with...\n');
fprintf('    --> %d terms in the expansion\n',cap_P);
fprintf('    --> highest overall polynomial degree: %d \n',p);
fprintf('    --> Total # of possible param combinations: %d\n',length(poly_roots)^N);
fprintf('    --> subsamped param combinations down to %d\n\n',N_subset);


%------------------------------------------------------------
% MULTIVARIABLE LEGENDRE POLYNOMIAL ORDERING
%------------------------------------------------------------
alphaMAT = create_Polynomial_Ordering(N,p);


%------------------------------------------------------------
% CREATE INFORMATION MATRIX
%------------------------------------------------------------
INFO_MAT = create_Information_Matrix(N,p,cap_P,param_combo_full,alphaMAT);

[matR,matC]=size(INFO_MAT);
fprintf('Information Matrix:\n');
fprintf('    --> Size: %dx%d\n',matR,matC);
fprintf('    --> Rank: %.8f\n\n',rank( INFO_MAT'*INFO_MAT, 1e-14 ) );


%------------------------------------------------------------
% SAVE ALL RELEVANT gPC INFO TO MATLAB '.mat' FILE
%------------------------------------------------------------
if save_Info
    save('gPC_Info.mat','N','p','cap_P','N_subset','param_combo','alphaMAT','INFO_MAT');
end


%---------------------------------------------------------------
% EVALUATE USER-SPECIFIED MODEL FOR PARTICULAR PARAMETER COMBOS
%---------------------------------------------------------------
Y = zeros(length(param_combo(:,1)),1);
for i=1:length(param_combo(:,1))
    Y(i,1) = USER_SPECIFIED_MODEL( param_combo(i,:) );
end



%------------------------------------------------------------
% SOLVE FOR COEFFICIENTS VIA LEAST-SQUARES WITH SVD
%          SVD:  A'* A = U Sigma V^*
%------------------------------------------------------------
[U,Sigma,V] = svd(INFO_MAT'* INFO_MAT);
%
% PSEUDO-INVERSE THROUGH SVD
sCoeffs = V*inv(Sigma)*ctranspose(U)*INFO_MAT'*Y;
%
% PSEUDO-INVERSE THROUGH STANDARD MATRIX INVERSION
%sCoeffs = inv(INFO_MAT'*INFO_MAT)*INFO_MAT'*Y;


%------------------------------------------------------------------------
% COMPUTE ERRORS BETWEEN gPC SURROGATE AND TRUE MODEL FUNCTION
%   -> how well does gPC surrogate recover training dataset
%   -> how well does gPC surrogate predict across independent test data
%------------------------------------------------------------------------
compute_Validation_and_Testing_Error(N,sCoeffs,alphaMAT,param_combo)


%------------------------------------------------------------------
% TEST EXPANSION ON 2D SUBSPACE:
%       -> plot response surface on 2D parameter subspace
%       -> qualitative/quantitative comparison on that subspace
%------------------------------------------------------------------
if test_Flag == 1
    test_gPC_Expansion(sCoeffs,alphaMAT)
end


%------------------------------------------------------------
% GET SOBOL INDICES
%------------------------------------------------------------
[sFirst,sTot,s2nd,s123] = compute_Sobol_Indices(sCoeffs,alphaMAT);


%------------------------------------------------------------
% PRINT SOME GENERAL INFORMATION TO SCREEN
%------------------------------------------------------------
fprintf('--------------------------------------\n');
fprintf('\nsX (1ST ORDER SOBOL) INDICES\n');
for n=1:N
    str = ['s' num2str(n)];
    fprintf('    --> %s: %.4f \n',str,sFirst(n));
end
%
fprintf('\nsXX (2ND ORDER) SOBOL INDICES\n');
for n1=1:N
    for n2=n1+1:N
    str= ['s' num2str(n1) num2str(n2)];
    fprintf('    --> %s: %.4f\n',str,s2nd(n1,n2));
    end
end
%
fprintf('\ns123 SOBOL INDEX\n');
fprintf('    --> x123: %.4f\n',s123(1));
%
fprintf('\nTOTAL-ORDER SOBOL INDICES\n');
for n=1:N
    str = ['sT' num2str(n)];
    fprintf('    --> %s: %.4f \n',str,sTot(n));
end
