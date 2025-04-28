%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Get Training and Testing Data for Neural Network
%         
%   Author: Nick Battista 
%   Date: March 23, 2023
%   Institution: TCNJ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TRAINING_DATA,TRAIN_OUTPUT,TESTING_DATA,TEST_OUTPUT,minZ,maxZ] = get_Training_and_Test_Data(numInputs,NTrain,NTest,epsError)

%------------------------------------------------------------------------
%
%             TRAINING DATA: Setting Up TRAINING SAMPLES
%
%------------------------------------------------------------------------

%-------------------------------
% SOBOL SEQUENCE ATTRIBUTES
%-------------------------------
skippy=0;
leapy=0;

%-----------------------------------------------
% SOBOL SEQUENCE:
%    -> get Sobol' sequence in [0,1]^N 
%    -> scale Sobol' sequence to [-1,1]^N
%-----------------------------------------------
sobol = sobolset(numInputs,'Skip',skippy,'Leap',leapy);
%
TRAINING_DATA = net(sobol,NTrain);     % data in [0,1]
%TRAINING_DATA = 2*net(sobol,NPts)-1; % data in [-1,1]
   

%------------------------------------------------------------------------
%
%             TESTING DATA: Set up TESTING SAMPLES
%
%------------------------------------------------------------------------

%----------------------------------------------------
% Find TEST Data Combinations
%----------------------------------------------------
TESTING_DATA = rand(NTest,numInputs);

%------------------------------------------------------------------------
%
%                     EVALUATE SPECIFIED FUNCTION 
%         (TO GET OUTPUT VALUES FOR TRAINING AND TESTING DATA)
%
%------------------------------------------------------------------------
flagScaleOutput = 0;
TRAIN_OUTPUT = Evaluate_Function(TRAINING_DATA,epsError,flagScaleOutput,0,0);
TEST_OUTPUT =  Evaluate_Function( TESTING_DATA,epsError,flagScaleOutput,0,0);

%-----------------------------------------
% SCALE OUTPUT DATA TO [0,1]
%-----------------------------------------
flagGetMinMax = 1;
if flagGetMinMax
    
    %-----------------------------------------
    % IDEA: y=mx+b
    %       y(minZ) = m(minZ) + b = 0;
    %       y(maxZ) = m(minZ) + b = 1; 
    %-----------------------------------------
    minZ = min([TRAIN_OUTPUT; TEST_OUTPUT]);
    maxZ = max([TRAIN_OUTPUT; TEST_OUTPUT]);
    %
    m = 1 / (maxZ-minZ);
    b = -m*minZ;
    %
    % Scale to [0,1]
    TRAIN_OUTPUT = m * TRAIN_OUTPUT+ b;
    TEST_OUTPUT = m *  TEST_OUTPUT+ b;
   
end

    