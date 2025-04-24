%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: compute Validation and Testing Error across full 3D parameter
%           space
%
%        -> Samples full space using Sobol' sequence
%
%       INPUTS:       N: # of uncertain parameters
%               sCoeffs: gPC expansion coefficients
%              alphaMAT: Legendre polynomial orderings
%           param_combo: sampled parameter combinations for training
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compute_Validation_and_Testing_Error(N,sCoeffs,alphaMAT,training_params)

    %-------------------------------
    % SOBOL SEQUENCE ATTRIBUTES
    %-------------------------------
    skippy=0;
    leapy=0;
    NPts = 1000;

    %-----------------------------------------------
    % SOBOL SEQUENCE:
    %    -> get Sobol' sequence in [0,1]^N 
    %    -> scale Sobol' sequence to [-1,1]^N
    %-----------------------------------------------
    p = sobolset(N,'Skip',skippy,'Leap',leapy);
    testing_params = 2*net(p,NPts)-1;
   
    
    %-----------------------------------------------
    % TRAINING DATA RECOVERY
    %-----------------------------------------------
    [gPC_Train,REAL_Train] = compute_gPC_and_TRUE_function_values(training_params,sCoeffs,alphaMAT);
    
    
    %-----------------------------------------------
    % TESTING DATA RECOVERY
    %-----------------------------------------------
    [gPC_Test,REAL_Test] = compute_gPC_and_TRUE_function_values(testing_params,sCoeffs,alphaMAT);

    %------------------------------------------------
    % COMPUTE ERROR:
    %    --> Resize data to column vector
    %    --> Compute percent error
    %------------------------------------------------
    errRel_Train = abs(gPC_Train-REAL_Train) ./ abs(REAL_Train) * 100;
    errRel_Test = abs(gPC_Test-REAL_Test) ./ abs(REAL_Test) * 100;
    
    
    fprintf('------------------------------------------\n\n');
    fprintf('gPC TRAINING Percent Error (Full 3D space):\n');
    fprintf('   --> AVG:    %.3f\n',mean( errRel_Train)  );
    fprintf('   --> MEDIAN: %.3f\n',median( errRel_Train)  );
    fprintf('   --> MIN:    %.3f\n',min( errRel_Train)  );
    fprintf('   --> MAX:    %.3f\n',max( errRel_Train)  );
    fprintf('   --> STDEV:  %.3f\n\n',std( errRel_Train)  );
    
    fprintf('------------------------------------------\n\n');
    fprintf('gPC TESTING Percent Error (Full 3D space):\n');
    fprintf('   --> AVG:    %.3f\n',mean( errRel_Test)  );
    fprintf('   --> MEDIAN: %.3f\n',median( errRel_Test)  );
    fprintf('   --> MIN:    %.3f\n',min( errRel_Test)  );
    fprintf('   --> MAX:    %.3f\n',max( errRel_Test)  );
    fprintf('   --> STDEV:  %.3f\n\n',std( errRel_Test)  );
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: compute the gPC predicted value and real function value 
%           at all parameter combinations in MAT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gPC_DATA,REAL_DATA] = compute_gPC_and_TRUE_function_values(MAT,sCoeffs,alphaMAT)

    % Allocate data storage
    gPC_DATA = zeros( length(MAT(:,1)), 1 );
    REAL_DATA = zeros( length(MAT(:,1)), 1 );

    % Loop over each parameter combination
    for n=1:length(MAT(:,1))
        
        % Evaluate Multidimensional Legendre Polynomial
        gPC_DATA(n,1) = evaluate_MultiDim_Legendre_Poly(sCoeffs, alphaMAT, MAT(n,:) ); 

        % TEST PROBLEM
        REAL_DATA(n,1) =  USER_SPECIFIED_MODEL( MAT(n,:) );
  
    end