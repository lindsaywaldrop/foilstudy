%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: trains, validates, and tests an Artificial Neural Network (ANN) 
%                        with *TWO* Hidden Layers and a bias in every layer!
%
% Author: Nick Battista
% Institution: TCNJ
% Date: 03/21/2022
%
% USER CAN VARY: 
%                [1] which function we're using a neural network to fit
%                           (in Evaluate_Function)
%                           
%                [2] amount of training data or test data used 
%
%                [3] HYPERPARAMETERS of the neural network, including:  
%                           - number of neurons in each hidden layer
%                           - learning rate (gradient descent stepsize)
%                           - number of epochs (full cycles of stochastic gradient
%                                   descent performed across all training data)
%                           - batch_size: how big each batch is in stochastic 
%                                   gradient descent
%
% SOME IDEAS: [1] The more training data you use, the longer it takes to
%                   'train' the neural network (same with more neurons)
%
%             [2] Using too many neurons can lead to 'overfitting', where
%                   the network recreates the training data very well but does
%                   not work well on test data
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function go_Go_ANN()

%---------------------------------------------------------------------------
% ADD PATH TO NN SOURCE SCRIPTS
%---------------------------------------------------------------------------
addpath('../src/');

%---------------------------------------------------------------------------
%                 GET TRAINING AND TESTING DATA:
%       --> Decide how you're sampling for each
%       --> Decide # of training samples
%       --> Chose a particular function in Evaluate_Function.m
%---------------------------------------------------------------------------
%
numInputs = 3;    % How many inputs  (should match size below)!
numOutputs = 1;   % How many outputs (should match size below)!
%
NTrain = 500;    % How many training samples
NTest = 500;     % How many testing samples
epsError = 0;    % How much epsilon percent error to add on average 
%
[TRAINING_DATA,TRAIN_OUTPUT,TESTING_DATA,TEST_OUTPUT,minZ,maxZ] = get_Training_and_Test_Data(numInputs,NTrain,NTest,epsError);

%
size(TRAINING_DATA)
size(TRAIN_OUTPUT)
%
size(TESTING_DATA)
size(TEST_OUTPUT)

%--------------------------------------------------------------------------
%       CONCATENATE TRAINING INPUTS AND LOG-TRANFORMED OUTPUTS
%--------------------------------------------------------------------------
%TRAINING_DATA = [TRAINING_DATA log10(TRAIN_OUTPUT)];
%TESTING_DATA = [TESTING_DATA log10(TEST_OUTPUT)];
%
TRAINING_DATA = [TRAINING_DATA (TRAIN_OUTPUT)]; % NO LOG-TRANSFORM OF OUTPUT
TESTING_DATA = [TESTING_DATA (TEST_OUTPUT)];    % NO LOG-TRANSFORM OF OUTPUT
%
fprintf('\n\n---------------------------------------------------\n');
fprintf('Train & Test Neural Network for different dimensions...\n');
fprintf('      --> all data in [0,1]\n');


%--------------------------------------------------------------------------
%        Save TRAINING and TEST Data (saves output in log form )
%--------------------------------------------------------------------------
strMAT = ['../src/matlab/nn/Training_and_Test_Data_eps' num2str(epsError) '.mat'];
save(strMAT,'TRAINING_DATA','TESTING_DATA');


%--------------------------------------------------------------------------
%           Transform the Training & Testing Data Inputs to [0,1]
%                          (if not already)
%--------------------------------------------------------------------------
%[TRAINING_DATA_SCALED,TESTING_DATA_SCALED,min_Z,max_Z] = scale_Training_and_Testing_Data_To_0_1(numOutputs,TRAINING_DATA,TESTING_DATA);
%
TRAINING_DATA_SCALED = TRAINING_DATA;
TESTING_DATA_SCALED = TESTING_DATA;


%--------------------------------------------------------------------------
%
%     -------- > > > > HYPERPARAMETER SELECTION < < < < --------
%
%--------------------------------------------------------------------------
num_hidden_layer_neurons = 200; % # of neurons per layer (size of hidden layer)
learning_rate = 0.55;           % learning rate ("lambda", Gradient Descent Step-size)
momentum = 0.0;                 % momentum ("inertia factor" from previous Gradient Descent Step)
maxEpochs = 1500;               % max. number of EPOCHS ( forward prop + back prop in training through ALL Training Data)
batch_size = 1;                 % # of samples per mini-batch in Stochastic Gradient Descent
adaptive_stepSize = 0;          % 0 for no, 1 for yes (Barzilai_Borwein step)
regularization_Flag = 0;        % 0 for none, 1 for L1, 2 for L2
lam_regularize = 8e-08;         % Regularization coefficient, lambda
%
% For Printing To Screen Purposes
print_Interval = 100;           % how often to print COST/ERROR info to screen during training
%
% Store hyperparameters into vector for input into TRAINING script
hyper_param_vec = [num_hidden_layer_neurons, learning_rate, momentum,...
                    maxEpochs, batch_size, print_Interval,...
                    adaptive_stepSize, numInputs,...
                    regularization_Flag, lam_regularize];




%--------------------------------------------------------------------------
%                   Train the ARTIFICIAL NEURAL NETWORK
%       
%       Returns: [1] All Weight Matrices
%                [2] All Bias Vectors
%                [3] Cost (error) vector
%--------------------------------------------------------------------------
fprintf('\n--> Starting ANN Training...\n');
[W1,W2,WEnd,b1,b2,bEnd,costVec,min_Cost] = Train_Artificial_Neural_Network(TRAINING_DATA_SCALED,TESTING_DATA_SCALED,hyper_param_vec);

%---------------------------------------------------
% Plot the COST/LOSS function vs. ITERATION NUMBER
%---------------------------------------------------
figure(5)
%
ms = 40;   % MarkerSize
lw=5;      % LineWidth
fs=22;     % FontSize
%
loglog(1:length(costVec),costVec,'.-','MarkerSize',ms,'LineWidth',lw); hold on;
xlabel('Epoch Number'); 
ylabel('Cost (loss)');
set(gca,'FontSize',fs);


%-----------------------------------------------------------
%
% SAVE INFORMATION (PRINT TO TXT FILES and .MAT FILES)
%
%-----------------------------------------------------------
currDir = pwd;
dirSave = '../src/matlab/nn/Trained_Neural_Network_Weights';
mkdir(dirSave);
cd(dirSave);
%
save('Trained_Neural_Network.mat','b1','b2','bEnd','W1','W2','WEnd','costVec','num_hidden_layer_neurons');
%
% print_Matrix_To_Txt_File(b1,'b1');
% print_Matrix_To_Txt_File(b1,'b2');
% print_Matrix_To_Txt_File(W1,'W1');
% print_Matrix_To_Txt_File(W1,'W2');
% print_Matrix_To_Txt_File(WEnd,'WEnd');
% print_Matrix_To_Txt_File(costVec(1:end)','costVec');
%
cd(currDir); % Go back to Original Directory


%-----------------------------------------------------------------------
%                           MODEL VALIDATION....
%
%            PERFORM FORWARD PROPAGATION ON *TRAINING* DATA
%               TO GET MODEL PREDICTED OUTPUT VALUES!
%                 (try to recovery training dataset)
%-----------------------------------------------------------------------
[zValid_Scaled,x2,x1] = forward_propagate(TRAINING_DATA_SCALED(:,1:numInputs)',W1,W2,WEnd,b1,b2,bEnd);
zValid_Scaled = zValid_Scaled';


%-----------------------------------------------------------------------
%                           MODEL TESTING....
%
%              PERFORM FORWARD PROPAGATION ON *TEST* DATA
%                 TO GET MODEL PREDICTED OUTPUT VALUES!
%             (to test model against independent test data)
%-----------------------------------------------------------------------
[zTest_Scaled,x2,x1] = forward_propagate(TESTING_DATA_SCALED(:,1:numInputs)',W1,W2,WEnd,b1,b2,bEnd);
zTest_Scaled = zTest_Scaled'; % make column vector
%
zValid = zValid_Scaled;
zTest = zTest_Scaled;

%-----------------------------------------------------------------------
%      UNDO LOGARITHMIC TRANSFORMATION (IF TRAINING WITH LOGARITHMS)!
%-----------------------------------------------------------------------
TRAINING_DATA(:,end) = TRAIN_OUTPUT;
TESTING_DATA(:,end) = TEST_OUTPUT;
%
%zValid = 10.^zValid_Scaled;
%zTest = 10.^zTest_Scaled;




%--------------------------------------------------------------------------
%
%         PRINT ERROR INFO (both TRAINING/VALIDATION & TEST DATA)
%
%--------------------------------------------------------------------------
%
% MODEL VALIDATION (how well did model do for TRAINING data)
print_Error_Information(zValid,TRAINING_DATA,'VALIDATION',numInputs);
%
% MODEL TESTING (how well did model do for TESTING data)
print_Error_Information(zTest,TESTING_DATA,'TESTING',numInputs);


%--------------------------------------------------------------------------
%                        MODEL VALIDATION PLOT: 
%    --> compare true values against Neural Network Predicted values!
%    --> qualitatively assess how well ANN recovers TRAINING DATA
%---------------------------------------------------------------------------
f10 = figure(10);
subplot(1,2,1)
plot(1:length(TRAINING_DATA(:,1)), TRAINING_DATA(:,numInputs+1),'k.','MarkerSize',ms+15); hold on;
plot(1:length(TRAINING_DATA(:,1)), zValid,'.','MarkerSize',ms,'Color',[0.2 0.675 0.75]); hold on;
leg=legend('Training Data','ANN Prediction (Validation)','Location','SouthEast');
xlabel('Data Point');
ylabel('Function Value: f(x1,x2,x3)');
axis([0 length(TRAINING_DATA(:,1))+1 0 1.05*max(TRAINING_DATA(:,numInputs+1))]);
title('ANN Model: Recover Training Dataset');
grid on;
grid minor;
set(gca,'FontSize',fs);
set(leg,'FontSize',20);

%--------------------------------------------------------------------------
%                           MODEL TESTING PLOT: 
%    --> compare true values against Neural Network Predicted values!
%    --> qualitatively assess how well ANN predicts independent TEST DATA
%--------------------------------------------------------------------------
subplot(1,2,2)
plot(1:length(TESTING_DATA(:,1)), TESTING_DATA(:,numInputs+1),'k.','MarkerSize',ms+15); hold on;
plot(1:length(TESTING_DATA(:,1)), zTest,'r.','MarkerSize',ms,'Color',[0.2 0.675 0.75]); hold on;
leg=legend('Testing Data','ANN Prediction (Testing)','Location','SouthEast');
xlabel('Data Point');
ylabel('Function Value: f(x1,x2,x3)');
axis([0 length(TESTING_DATA(:,1))+1 0 1.05*max(TESTING_DATA(:,numInputs+1))]);
title('ANN Model: Validation/Testing');
grid on;
grid minor;
set(gca,'FontSize',fs);
set(leg,'FontSize',20);
%
f10.Position = [100 100 1800 450];


%-----------------------------------------------------------------------
%
%               Test 2-D Slice Through overall domain!
%              
%-----------------------------------------------------------------------
x1_Vec = linspace(0,1,100); 
x2_Vec = linspace(0,1,100);
[X1,X2]=meshgrid(x1_Vec,x2_Vec);
%
% Allocate storage for data
zPRED = zeros(size(X1));
zREAL = zPRED;
%
x3 = 0.5;  % CONSTANT VALUE OF THIRD PARAMETER
%
% Loop thru all values of x1 in x1_Vec
for n1=1:length(x1_Vec)
    
    % Get Specific x1 input value
    x1 = x1_Vec(n1);
    
    % Loop thru all entries of x2 in x2_Vec
    for n2=1:length(x2_Vec)
        
        % Get Specific x2 input value
        x2=x2_Vec(n2);
        
        % Get particular (x1,x2,x3) combination to pass through NN
        InputVec = [x1 x2 x3]';
        
        % Get Predicted Output 
        [zOUT,~,~] = forward_propagate(InputVec,W1,W2,WEnd,b1,b2,bEnd);
        
        % Un-log-transform output (if NN was trained on logarithmic data)
        %zPRED(n2,n1) = 10.^(zOUT);
        
        % Predicted Data Value
        zPRED(n2,n1) = zOUT;
        
        % Get True Data Values
        zREAL(n2,n1) = Evaluate_Function(InputVec',eps,1,minZ,maxZ);
        
    end 
end

%-----------------------------------------------------------------------
%                           PLOT 2-D SLICE!
%-----------------------------------------------------------------------
f55=figure(55);
%
subplot(1,3,1)
surf(X1,X2,zREAL)
xlabel('x1');
ylabel('x2');
zlabel('f(x1,x2,x3) [TRUE]');
set(gca,'FontSize',fs);
zlim([0 1]);
%
subplot(1,3,2)
surf(X1,X2,zPRED)
xlabel('x1');
ylabel('x2');
zlabel('f(x1,x2,x3) [ANN PREDICTION]');
set(gca,'FontSize',fs);
zlim([0 1]);
%
subplot(1,3,3)
surf(X1,X2,abs(zREAL-zPRED)./zREAL);
xlabel('x1');
ylabel('x2');
zlabel('Relative Error: |Real-ANN|/Real');
set(gca,'FontSize',fs);

f55.Position = [50 500 1600 450];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Scale the TRAINING AND TESTING DATA for all quantities to be [0,1]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TRAINING_DATA_SCALED,TESTING_DATA_SCALED,min_Z,max_Z] = scale_Training_and_Testing_Data_To_0_1(numOutputs,TRAINING_DATA,TESTING_DATA)

%-----------------------------------------------------
% Allocated Storage
%-----------------------------------------------------
%TRAINING_DATA_SCALED = zeros( size( TRAINING_DATA ) );
%TESTING_DATA_SCALED = zeros( size( TESTING_DATA ) );
%
TRAINING_DATA_SCALED = TRAINING_DATA;
TESTING_DATA_SCALED = TESTING_DATA;

%-----------------------------------------------------
%
%      >>>>> Desired 5D Parameter Space <<<<<
%
%-----------------------------------------------------
%
% REYNOLDS NUMBER
p1_min = 25;
p1_max = 275;

% CONTRACTION PERIOD
p2_min = 0.07;
p2_max = 0.43;

% EXPANSION PERIOD
p3_min = 0.38;
p3_max = 0.82;

% REST-2 PERIOD
p4_min = 0.0;
p4_max = 0.42;

% BENDING STIFFNESS (NOTE THE LOG10)
p5_min = log10(2e4);
p5_max = log10(8e5);

% STORE VALUES IN VECTOR FOR AUTOMATED SCALING
SCALE_MIN = [p1_min p2_min p3_min p4_min p5_min];
SCALE_MAX = [p1_max p2_max p3_max p4_max p5_max];


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
% FUNCTION prints matrix to file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Matrix_To_Txt_File(a,strName)

nameTxt = [strName '.txt'];

fid = fopen(nameTxt, 'wt'); % Open for writing
for i=1:size(a,1)
   fprintf(fid, '%.12f ', a(i,:));
   fprintf(fid, '\n');
end
fclose(fid);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: performs/prints Error Information for Model Validation & Model
%           Testing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Error_Information(zData,DATA,strKind,numInputs)


% FIRST COMPUTE ERRORS
errABS = abs( zData-DATA(:,numInputs+1) );       % Calculate Absolute Error
errREL = ( errABS ./ DATA(:,numInputs+1) ) * 100; % Calculate Relative Error Percentage
%
inds = isinf(errREL); % finds inds if divide by zero occurs
errREL(inds) = 555;   % if divide by zero happens, revert to 100% error.
%                     %        NOTE: will happen if yData = 0
%
% THEN PRINT INFO TO SCREEN
%
stringString = ['*** ' strKind ' ERROR INFORMATION ***\n'];
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
%



