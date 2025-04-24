function Prediction_Errors()

addpath('../');

%------------------------------------------
% LOADS: neural network weights!
%       --> W1, W2, b1, b2, WEnd
%------------------------------------------
load('Trained_Neural_Network_Weights_Eps_0pt10/Trained_Neural_Network.mat');

numInputs = 7;
numOutputs = 1;

%---------------------------------------------------------------------------
% GET TRAINING AND TESTING DATA: 
%           NOTE --> gives output already in log10 form
%---------------------------------------------------------------------------
load('Training_and_Test_Data_0pt10.mat');
size(TRAINING_DATA)
size(TESTING_DATA)
%
%
fprintf('Using PERFECT Test Data to Compare\n');
epsError = 0.0;
TEST_OUTPUT = Evaluate_Model(TESTING_DATA,epsError);
TESTING_DATA(:,end) = log10(TEST_OUTPUT);


%-------------------------------------------------------
% Scale the Training and Testing Data to [0,1]
%-------------------------------------------------------
%[TRAINING_DATA_SCALED,TESTING_DATA_SCALED,min_Z,max_Z] = scale_Training_and_Testing_Data_To_0_1(numOutputs,TRAINING_DATA,TESTING_DATA);
%
TRAINING_DATA_SCALED = TRAINING_DATA;
TESTING_DATA_SCALED = TESTING_DATA;

%--------------------------------------------------
%             MODEL VALIDATION....
%
% PERFORM FORWARD PROPAGATION ON *TRAINING* DATA
%    TO GET MODEL PREDICTED OUTPUT VALUES!
%---------------------------------------------------
[zValid_Scaled,x2,x1] = forward_propagate(TRAINING_DATA_SCALED(:,1:numInputs)',W1,W2,WEnd,b1,b2);
zValid_Scaled = zValid_Scaled';



%----------------------------------------------
%              MODEL TESTING....
%
% PERFORM FORWARD PROPAGATION ON *TEST* DATA
%    TO GET MODEL PREDICTED OUTPUT VALUES!
%-----------------------------------------------
[zTest_Scaled,x2,x1] = forward_propagate(TESTING_DATA_SCALED(:,1:numInputs)',W1,W2,WEnd,b1,b2);
zTest_Scaled = zTest_Scaled'; % make column vector


%-----------------------------------------------
% Transform appropriately!
%-----------------------------------------------
%
% IF NOT USING LOGS IN ERROR CALCULATIONS
TRAINING_DATA(:,numInputs+1:end) = 10.^( TRAINING_DATA(:,numInputs+1:end) );
TESTING_DATA(:,numInputs+1:end) = 10.^( TESTING_DATA(:,numInputs+1:end) );
%
zValid_Scaled = 10.^( zValid_Scaled );
zTest_Scaled = 10.^( zTest_Scaled );

%--------------------------------------------

zValid = zValid_Scaled;
zTest = zTest_Scaled;

% TRAINING_DATA(:,numInputs+1:end) = 10.^TRAINING_DATA(:,numInputs+1:end);
% TESTING_DATA(:,numInputs+1:end) = 10.^TESTING_DATA(:,numInputs+1:end);

% mean( abs(zValid(:,1)-TRAINING_DATA(:,5))./TRAINING_DATA(:,5) )
% mean( abs(zValid(:,2)-TRAINING_DATA(:,6))./TRAINING_DATA(:,6) )
% mean( abs(zValid(:,3)-TRAINING_DATA(:,7))./TRAINING_DATA(:,7) )
% 
% asga

%--------------------------------------------------------------------------
% SCALE output data back to REAL output range, i.e., 
%                                           from [0,1]->[min_Z,max_Z]
%--------------------------------------------------------------------------
%[zValid,zTest] = scale_Back_To_Real_Output_Range(min_Z,max_Z,zValid_Scaled,zTest_Scaled);


strVec = {'Output'};
for j=1:numOutputs

    %-----------------------------------------------------------------------------
    % MODEL VALIDATION PLOT: jth PERFORMANCE METRIC-j
    %       --> compare true values against Neural Network Predicted values!
    %       --> "see how well ANN did on TRAINING data"
    %-----------------------------------------------------------------------------
    %
    ms = 30;   % MarkerSize
    lw=5;      % LineWidth
    fs=22;     % FontSize
    %
    fJ = figure(j);
    subplot(1,2,1)
    plot(1:length(TRAINING_DATA(:,1)), TRAINING_DATA(:,numInputs+j),'k.','MarkerSize',ms+15); hold on;
    plot(1:length(TRAINING_DATA(:,1)), zValid(:,j),'r.','MarkerSize',ms); hold on;
    leg=legend('Training Data','ANN Prediction (Validation)');
    xlabel('Data Point');
    ylabel('z Value');
    %axis([0 length(TRAINING_DATA(:,1))+1 0 1.05*max(TRAINING_DATA(:,numInputs+j))]);
    strTitle = [strVec{j} ' Validation'];
    title(strTitle);
    set(gca,'FontSize',fs);
    set(leg,'FontSize',18);

    %-----------------------------------------------------------------------------
    % MODEL TESTING PLOT: jth PERFORMANCE METRIC
    %       --> compare true values against Neural Network Predicted values!
    %       --> "see how well ANN did on TESTING data"
    %-----------------------------------------------------------------------------
    subplot(1,2,2)
    plot(1:length(TESTING_DATA(:,1)), TESTING_DATA(:,numInputs+j),'k.','MarkerSize',ms+15); hold on;
    plot(1:length(TESTING_DATA(:,1)), zTest(:,j),'r.','MarkerSize',ms); hold on;xlabel('Speed');
    leg=legend('Testing Data','ANN Prediction (Testing)');
    xlabel('Data Point');
    ylabel('z Value');
    %axis([0 length(TESTING_DATA(:,1))+1 0 1.05*max(TESTING_DATA(:,numInputs+1))]);
    strTitle = [strVec{j} ' Testing'];
    title(strTitle);
    set(gca,'FontSize',fs);
    set(leg,'FontSize',18);
    %
    fJ.Position = [100 100 1800 450];

end

%------------------------------------------------------------
%
% PRINT ERROR INFO (both TRAINING/VALIDATION & TEST DATA)
%
%------------------------------------------------------------
%
% MODEL VALIDATION (how well did model do for TRAINING data)
print_Error_Information(zValid,TRAINING_DATA,'VALIDATION',numInputs,numOutputs);
%
% MODEL TESTING (how well did model do for TESTING data)
print_Error_Information(zTest,TESTING_DATA,'TESTING',numInputs,numOutputs);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Scale the TRAINING AND TESTING DATA for all quantities to be [0,1]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TRAINING_DATA_SCALED,TESTING_DATA_SCALED,min_Z,max_Z] = scale_Training_and_Testing_Data_To_0_1(numOutputs,TRAINING_DATA,TESTING_DATA)

%-----------------------------------------------------
% Allocated Storage
%-----------------------------------------------------
% TRAINING_DATA_SCALED = zeros( size( TRAINING_DATA ) );
% TESTING_DATA_SCALED = zeros( size( TESTING_DATA ) );
%
TRAINING_DATA_SCALED = TRAINING_DATA;
TESTING_DATA_SCALED = TESTING_DATA;

%-----------------------------------------------------
%
%      >>>>> Desired 2D Parameter Space <<<<<
%
%-----------------------------------------------------
%
% LEG FRAC
p1_min = 0.55;
p1_max = 0.95;

% Leg NUM
p2_min = 10;
p2_max = 30;

% % TENTACLE STIFFNESS
% p3_min = log10(3e7);
% p3_max = log10(3e8);
% 
% % CONTRACTION AMOUNT (based on angle)
% p4_min = -17*pi/180;
% p4_max =  35*pi/180;


% STORE VALUES IN VECTOR FOR AUTOMATED SCALING
SCALE_MIN = [p1_min p2_min];% p3_min p4_min];
SCALE_MAX = [p1_max p2_max];% p3_max p4_max];


%---------------------------------------------------
% Scale INPUT values from [min,max]->[0,1]
%---------------------------------------------------
for i=1:length(TRAINING_DATA(1,:))-numOutputs

    %-----------------------------------------
    % Scale to [0,1]: [1]  min_Z*m + b = 0
    %                 [2]  max_Z*m + b = 1
    %-----------------------------------------
    min_Z = SCALE_MIN(i);
    max_Z = SCALE_MAX(i);
    %
    m = 1/(max_Z-min_Z);
    b = 1-max_Z*m;
    %
    TRAINING_DATA_SCALED(:,i) = m * TRAINING_DATA(:,i) + b;
    TESTING_DATA_SCALED(:,i) =  m * TESTING_DATA(:,i) + b;

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Scale Neural Network output from [0,1] to *REAL* output range
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [zValid,zTest] = scale_Back_To_Real_Output_Range(min_Z,max_Z,zValid_Scaled,zTest_Scaled)
%    
%IDEA: y=mx+b
%      y(0) =   0*m + b = min_Z;
%      y(1)  =   1m + b = max_Z; 
%
% Linear Scaling Variables
b = min_Z;
m = max_Z-min_Z;
%
%
% Scale Back to Original Data Levels (for predicted validation & test data)
zValid = m * zValid_Scaled + b;
zTest = m * zTest_Scaled + b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: performs/prints Error Information for Model Validation & Model
%           Testing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Error_Information(zData,DATA,strKind,numInputs,numOutputs)

strVec = {'Output'};

for n = 1:numOutputs
    
    nn = numInputs+n;
    
    % FIRST COMPUTE ERRORS
    errABS = abs( zData(:,n) - DATA(:,nn) );           % Calculate Absolute Error
    errREL = ( errABS ./ DATA(:,nn) ) * 100; % Calculate Relative Error Percentage
    %
    inds = isinf(errREL); % finds inds if divide by zero occurs
    errREL(inds) = 555;   % if divide by zero happens, revert to 100% error.
    %                     %        NOTE: will happen if yData = 0
    %
    % THEN PRINT INFO TO SCREEN
    %
    stringString = ['*** ' strVec{n} ' ' strKind ' ERROR INFORMATION ***\n'];
    %
    fprintf('\n\n***---------------------------***\n');
    fprintf(stringString);
    fprintf('***---------------------------***\n\n');
    %
    fprintf(' --> ABSOLUTE ERROR  <-- \n');
    fprintf('Max ABS Error: %.4f\n', max( errABS(:,1) ) );
    fprintf('Min ABS Error: %.4f\n', min( errABS(:,1) ) );
    fprintf('Avg ABS Error: %.4f\n', mean( errABS(:,1) ) );
    fprintf('\n');
    %
    fprintf(' --> RELATIVE PERCENT ERROR  <-- \n');
    fprintf('Max Percent Error: %.4f\n', max( errREL(:,1) ) );
    fprintf('Min Percent Error: %.4f\n', min( errREL(:,1) ) );
    fprintf('Avg Percent Error: %.4f\n', mean( errREL(:,1) ) );
    fprintf('\n');

    [maEx,ind]=max( errREL(:,1) )
    
    clear errABS
    clear errREL
    
end
