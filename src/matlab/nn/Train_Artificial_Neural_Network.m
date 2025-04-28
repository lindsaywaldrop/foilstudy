%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: TRAINS an artificial neural network with 2 hidden layers
%
%  Hyperparameters: (1) num_hidden_layer_neurons
%   
%                   (2) maxIter: times to iterate gradient descent
%                                               
%                   (3) regularize_Flag: choice of regularization
%                                                     
%                   (4) lam_regularize: penalty-term for regularization
%                                                    
%                   (5) learning rates: lambda_1, lambda_2,...
%                       NOTE: code tries to do adaptive learning rate but 
%                             can define a minimum learning  rate (minLAM)
%                                        
%                   (6) batch_size: amount of data to use in each mini
%                                   batch of stochastic gradient descent
%                    
%                   (7) momentum: **NOT INCORPORATED HERE**
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W1_Save,W2_Save,WEnd_Save,b1_Save,b2_Save,bEnd_Save,...
                    costVec,costVecTrain,costVecTest,min_Cost] = Train_Artificial_Neural_Network(TRAINING_DATA,TESTING_DATA_SCALED,hyper_param_vec)

%------------------------------------------------------
% FLAG: use previous NN weights
%------------------------------------------------------
flag_Use_Stored_Weights = 0;
%
if flag_Use_Stored_Weights
    fprintf('\n   --> Loading previously stored weights...\n');
    fprintf('             (Trying to make them more accurate...)\n\n');
else
    fprintf('\n   --> NOT starting w/ previously stored weights...\n\n');
end

%------------------------------------------------------
% hyper_param_vec(1): num_hidden_layer_neurons
% hyper_param_vec(2): learning_rate
% hyper_param_vec(3): momentum
% hyper_param_vec(4): number of epochs
% hyper_param_vec(5): batch_size
% hyper_param_vec(6): print_Interval
% hyper_param_vec(7): flag for adaptive step size
% hyper_param_vec(8): # of input parameters
% hyper_param_vec(9): flag for regularization (0-none, 1-L1, 2-L2)
% hyper_param_vec(10): regularization value
%------------------------------------------------------

%----------------------------------
% Redefine Data for Consistency
%----------------------------------
%
x0 = TRAINING_DATA(:,1:hyper_param_vec(8))'; % should be of size [ - , numData]; - <--- number of pts per single data
z0 = TRAINING_DATA(:,hyper_param_vec(8)+1:end)'; % should be of size [ -, numData ];
%
num_input_neurons = hyper_param_vec(8);
train_data_size = length( x0(1,:) );
%
x0_TEST = TESTING_DATA_SCALED(:,1:hyper_param_vec(8))';
z0_TEST = TESTING_DATA_SCALED(:,hyper_param_vec(8)+1:end)';



%----------------------------------------------
% Rule of thumb for # OF HIDDEN LAYER NEURONS
%               (Hidden Layer Size)
%      (https://stats.stackexchange.com/questions/181/how-to-choose-the-number-of-hidden-layers-and-nodes-in-a-feedforward-neural-netw )
%           %alpha = 10;                 % Hyperparameter between [2,10]
%           %num_hid_layer_neurons = ceil( numTrain / ( alpha *( Ninput + Noutput ) ) )
%----------------------------------------------
Ninput = num_input_neurons;  % Amount of INPUT neurons
Noutput = length(z0(:,1));   % Amount of OUTPUT neurons
numTrain = train_data_size;  % Amount of Training Data
%
num_hid_layer_neurons = hyper_param_vec(1);

%
%---------------------------------------
% Get Number of DATA OUTPUTs
%---------------------------------------
num_Z = numel( z0 );         % # of data elements in Z (training set)
num_ZTest = numel(z0_TEST);  % # of data elements in Z (test set)

%---------------------------------------------------------------
%---------------------------------------------------------------
% Initialize Weights and Biases
%---------------------------------------------------------------
%---------------------------------------------------------------
if flag_Use_Stored_Weights 
    % Will load {W1,W2,WEnd} +  {b1,b2} + {num_hid_layer_neurons}
    load('Trained_Neural_Network_Weights/Trained_Neural_Network.mat');
    clear costVec;
else
    
    %----------------------------------------
    % INITIALIZE WEIGHTS/BIASES RANDOMLY
    %----------------------------------------
    coeff = 1 / sqrt( numTrain ); % for randomizing initial values
    W1 = coeff*( 2*rand( num_hid_layer_neurons, num_input_neurons)-1 );
    W2 = coeff*( 2*rand( num_hid_layer_neurons, num_hid_layer_neurons)-1 );
    WEnd = coeff*( 2*rand( Noutput, num_hid_layer_neurons )-1 );
    b1 = coeff*( 2*rand( num_hid_layer_neurons, 1 )-1 );   
    b2 = coeff*( 2*rand( num_hid_layer_neurons, 1 )-1 );
    bEnd = coeff*( 2*rand( Noutput, 1 )-1 );

end

%-------------------------------------------------------------
% INITIALIZE GRADIENTS OF WEIGHT/BIAS MATRICES ---------------
%-------------------------------------------------------------
W1_p = zeros( num_hid_layer_neurons, num_input_neurons );
dJ_dW1_p = zeros( num_hid_layer_neurons, num_input_neurons );
%
W2_p = zeros( num_hid_layer_neurons, num_hid_layer_neurons );
dJ_dW2_p = zeros( num_hid_layer_neurons, num_hid_layer_neurons );
%
WEnd_p = zeros( Noutput, num_hid_layer_neurons );
dJ_dWEnd_p = zeros( Noutput, num_hid_layer_neurons );
%
b1_p = zeros( num_hid_layer_neurons, 1 );
dJ_db1_p = zeros( num_hid_layer_neurons, 1 );
%
b2_p = zeros( num_hid_layer_neurons, 1 );
dJ_db2_p = zeros( num_hid_layer_neurons, 1 );
%
bEnd_p = zeros( num_hid_layer_neurons, 1 );
dJ_dbEnd_p = zeros( num_hid_layer_neurons, 1 );

%---------------------------------------------
% Initialize Learning Rates and Momentum
%---------------------------------------------
% Weight Learning Rates
lambda_1 = hyper_param_vec(2);     % initial learning rate
lambda_2 = hyper_param_vec(2);     % initial learning rate
lambda_End = hyper_param_vec(2);   % initial learning rate
% Bias Learning Rates
lambda_b1 = hyper_param_vec(2);    % initial learning rate
lambda_b2 = hyper_param_vec(2);    % initial learning rate
lambda_bEnd = hyper_param_vec(2);  % initial learning rate
% Minimum Lambda 
minLAM = hyper_param_vec(2);       % MINIMAL learning rate
%
alpha = hyper_param_vec(3);        % Momentum (*NOT USED*)
%
flag_Adaptive_StepSize = hyper_param_vec(7); % Adapative Step-Size w/ Barzilai_Borwein
lamCT = 0;
lamVec = [0.01 0.005 0.0025 0.001 0.0005 0.0002 0.0001 0.00005 0.00002];

%------------------------------------------------
% BATCH SIZE FOR MINIBATCH GRADIENT DESCENT
%           (vein of Stochastic Grad. Descent)
%------------------------------------------------
batch_size_Save = hyper_param_vec(5);


%---------------------------------------
% Initialize Regularization Variables
%---------------------------------------
lam_regularize = hyper_param_vec(10);
regFLAG = hyper_param_vec(9);
if regFLAG == 0
    regularize_Flag = 'none'; 
elseif regFLAG == 1
    regularize_Flag = 'L1';
elseif regFLAG == 2
    regularize_Flag = 'L2';
else
    regularize_Flag = 'none';
end
                        

                          
%-------------------------------------------------------------
%-------------------------------------------------------------
% START TRAINING!
%-------------------------------------------------------------
%-------------------------------------------------------------
numEpochs = hyper_param_vec(4);       % max # of EPOCHS to TRAIN model
costVec = zeros( 1, numEpochs );      % initalize vector storage for cost function
costVecTest = zeros( 1, numEpochs );  % initalize vector storage for cost function
costVecTrain = zeros( 1, numEpochs ); % initalize vector storage for cost function
%
for epochIter=1:numEpochs
    

    %-----------------------------------------------------------
    % PSEUDO-ADAPTIVE STEP-SIZE
    %-----------------------------------------------------------  
    if mod(epochIter,251)==0 && epochIter>=200 && flag_Adaptive_StepSize == 1
        lamCT = lamCT+1;

        % Scale learning rate
        minLAM = 0.90*minLAM;

        % Use stored learning rate values
        %minLAM = lamVec(lamCT);

        % Give computer time to rest
        %pause(1);

    % Just give computer time to rest    
    %elseif mod(epochIter,201)==0 && flag_Adaptive_StepSize == 1
    %   pause(1);
    end
    %
    lambda_1 = minLAM;
    lambda_2 = minLAM;
    lambda_End = minLAM;
    lambda_b1 = minLAM;
    lambda_b2 = minLAM;
    lambda_bEnd = minLAM;    
    
    %-----------------------------------------------------------
    % RANDOMLY SHUFFLE INDICES for MINI-BATCH RANDOM SAMPLING
    %       and RESET EPOCH PARAMETERS
    %-----------------------------------------------------------
    indsRandom = randperm(length(1:numTrain));  % Randomly shuffle training data indices for SGD
    costSum = 0;                                % Reset SINGLE Epoch cost to 0
    batch_size = batch_size_Save;               % Reset to original batch size
    
    %----------------------------------------------------------------
    % Iteration Number Inside SINGLE EPOCH
    %----------------------------------------------------------------
    for iter=1:floor(numTrain/batch_size_Save)


        %----------------------------------------------------------------
        % RANDOMLY SHUFFLE INDICES for MINI-BATCH RANDOM SAMPLING
        %----------------------------------------------------------------
        if iter ~= floor(numTrain/batch_size)
            inds = indsRandom( 1+(iter-1)*batch_size:batch_size*iter);
        else
            inds = indsRandom( 1+(iter-1)*batch_size:end);
            batch_size = length(inds);
        end


        %----------------------------------------------------------------
        % Forward Propagation
        %----------------------------------------------------------------
        [zHat,x2,x1] = forward_propagate(x0(:,inds),W1,W2,WEnd,b1,b2,bEnd);


        %----------------------------------------------------------------
        % Compute Cost Function (w/ or w/o regularization)
        %----------------------------------------------------------------
        %
        % Orig. Vectorized Mean-Squared Error (MSE) Cost ("loss") Function --------
        %
        J = 0.5*( z0(:,inds) - zHat ).^2;
        %
        % Regularize the Cost Function --------------------------------------------
        %
        if strcmp(regularize_Flag,'L2')
            cost =  ( sum(sum(J) ) + 0.5*lam_regularize*( sum( sum( W1.^2 ) )  + sum( sum( WEnd.^2 ) ) + sum( sum( b1.^2 ) )    ) );
        elseif strcmp(regularize_Flag,'L1') 
            cost = ( sum(sum(J) ) + 0.5*lam_regularize*( sum( sum( abs(W1) ) ) + sum( sum( abs( WEnd ) ) ) + sum( sum( abs(b1) ) )  ) ); 
        else
            cost = ( sum( sum(J) ) ); 
        end    
        %
        costSum = costSum + cost; % cumulative cost across epoch

        
        %----------------------------------------------------------------
        % Compute Delta Matrices (Fall 2022 w/ Eileen Yizzi)
        %       --> uses f_out = x; so f'_{out}=1
        %----------------------------------------------------------------
        delta_END = -( z0(:,inds) - zHat );
        delta_2 = (WEnd'*delta_END) .* act_function_PRIME( W2*x1 + b2 );
        delta_1 =     (W2'*delta_2) .* act_function_PRIME( W1*x0(:,inds) + b1 );
        
        
        %----------------------------------------------------------------
        % Compute REGULARIZATION Gradients
        %----------------------------------------------------------------
        if strcmp( regularize_Flag , 'L2' )
            regularize_grad_W1 = lam_regularize*( W1 );
            regularize_grad_W2 = lam_regularize*( W2 );
            regularize_grad_WEnd = lam_regularize*( WEnd );
            regularize_grad_b1 = lam_regularize*( b1 );
            regularize_grad_b2 = lam_regularize*( b2 );
            regularize_grad_bEnd = lam_regularize*( bEnd );            
        elseif strcmp( regularize_Flag , 'L1' )
            regularize_grad_W1 = lam_regularize*( sign(W1) );
            regularize_grad_W2 = lam_regularize*( sign(W2) );
            regularize_grad_WEnd = lam_regularize*( sign(WEnd) );
            regularize_grad_b1 = lam_regularize*( sign(b1) );
            regularize_grad_b2 = lam_regularize*( sign(b2) );
            regularize_grad_bEnd = lam_regularize*( sign(bEnd) );            
        else
            regularize_grad_W1 = 0;
            regularize_grad_W2 = 0;
            regularize_grad_WEnd = 0;
            regularize_grad_b1 = 0;
            regularize_grad_b2 = 0;
            regularize_grad_bEnd = 0;            
        end


        
        %----------------------------------------------------------------
        % COMPUTE GRADIENTS!  [>>> w/ Eileen Yizzi (Fall 2022) <<<]
        %
        %           Note: 1/num_Z added from COST function coeff.
        %                       (not included in partial deriv. derivation)
        %----------------------------------------------------------------
        %
        dJ_dWEnd = 1/batch_size * ( ( delta_END * x2' ) + regularize_grad_WEnd  );           
        %
        dJ_dW2 =   1/batch_size * ( delta_2 * x1'  + regularize_grad_W2  );
        %
        dJ_dW1 =   1/batch_size * ( delta_1 * x0(:,inds)'  + regularize_grad_W1  );
        %
        dJ_dbEnd = 1/batch_size * ( delta_END * ones(batch_size,1)  + regularize_grad_bEnd );        
        %
        dJ_db2 =   1/batch_size * ( delta_2 * ones(batch_size,1)  + regularize_grad_b2 );
        %
        dJ_db1 =   1/batch_size * ( delta_1 * ones(batch_size,1)  + regularize_grad_b1 );

 

        %------------------------------------------------------------------
        % Perform Gradient Descent (with momentum <<-- not incorporated)
        %------------------------------------------------------------------
        %
        % WEIGHTS!!! -------------------------------
        %
        W1n = W1 - lambda_1 * ( dJ_dW1 ); % - alpha*lambda_1*dJ_dW1_p;
        %
        W2n = W2 - lambda_2 * ( dJ_dW2 ); % - alpha*lambda_2*dJ_dW2_p;
        %
        WEndn = WEnd - lambda_End * ( dJ_dWEnd ); % - alpha*lambda_End*dJ_dWEnd_p;
        %
        % BIASES!!! --------------------------------
        %
        b1n = b1 - lambda_b1 * ( dJ_db1 ); % - alpha * lambda_b1 * dJ_db1_p;
        %
        b2n = b2 - lambda_b2 * ( dJ_db2 ); % - alpha * lambda_b2 * dJ_db2_p;
        %
        bEndn = bEnd - lambda_bEnd * ( dJ_dbEnd ); % - alpha * lambda_b2 * dJ_db2_p;

       
        %----------------------------------------------------------------
        % Use adaptive step-sizes (if flag = 1)
        %----------------------------------------------------------------
        if flag_Adaptive_StepSize
            lambda_1 = get_Barzilai_Borwein_Step( W1, W1_p, dJ_dW1, dJ_dW1_p );
            lambda_2 = get_Barzilai_Borwein_Step( W2, W2_p, dJ_dW2, dJ_dW2_p );
            lambda_End = get_Barzilai_Borwein_Step( WEnd, WEnd_p, dJ_dWEnd, dJ_dWEnd_p );
            %
            lambda_b1 = get_Barzilai_Borwein_Step( b1, b1_p, dJ_db1, dJ_db1_p );
            lambda_b2 = get_Barzilai_Borwein_Step( b2, b2_p, dJ_db2, dJ_db2_p );
            lambda_bEnd = get_Barzilai_Borwein_Step( bEnd, bEnd_p, dJ_dbEnd, dJ_dbEnd_p );            
            %
            if ( ( lambda_1 < minLAM ) || ( isnan( lambda_1 ) ) )
               lambda_1 = minLAM;
            end
            if ( ( lambda_2 < minLAM ) || ( isnan( lambda_2 ) ) )
               lambda_2 = minLAM;
            end
            if ( ( lambda_End < minLAM ) || ( isnan( lambda_End ) ) )
               lambda_End = minLAM;
            end
            if ( ( lambda_b1 < minLAM ) || ( isnan( lambda_b1 ) ) )
               lambda_b1 = minLAM;
            end
            if ( ( lambda_b2 < minLAM ) || ( isnan( lambda_b2 ) ) )
               lambda_b2 = minLAM;
            end
            if ( ( lambda_bEnd < minLAM ) || ( isnan( lambda_bEnd ) ) )
               lambda_bEnd = minLAM;
            end            
     
            %---------------------------------------------------------------
            % Save Previous Iteration data for Adaptive Gradient Step-Sizes
            %---------------------------------------------------------------
            W1_p = W1;
            W2_p = W2;
            WEnd_p = WEnd;
            %
            dJ_dW1_p = dJ_dW1;
            dJ_dW2_p = dJ_dW2;
            dJ_dWEnd_p = dJ_dWEnd;
            %
            b1_p = b1;
            b2_p = b2;
            bEnd_p = bEnd;            
            %
            dJ_db1_p = dJ_db1;
            dJ_db2_p = dJ_db2;
            dJ_dbEnd_p = dJ_dbEnd;            
            
        end

        %----------------------------------------------------------------
        % Increment to next
        %----------------------------------------------------------------
        W1 = W1n;
        W2 = W2n;
        WEnd = WEndn;
        %
        b1 = b1n;
        b2 = b2n;
        bEnd = bEndn;

        

    end
    
   
    %----------------------------------------------------------------
    %            Save COST for SINGLE Epoch (TRAINING)
    %----------------------------------------------------------------
    costVec(epochIter) = costSum;    % Store Cost For SINGLE EPOCH

    
    %-----------------------------------------------------------------
    %           Save COST for SINGLE Epoch (TESTING DATA)
    %-----------------------------------------------------------------
    [zHatTrain,~,~] = forward_propagate(x0,W1,W2,WEnd,b1,b2,bEnd);
    %
    J_Train = sum( sum( 0.5*( z0 - zHatTrain ).^2 ) );  
    %
    costVecTrain(epochIter) = J_Train / num_Z; % Store(TESTING) Cost For SINGLE EPOCH


    %-----------------------------------------------------------------
    %           Save COST for SINGLE Epoch (TESTING DATA)
    %-----------------------------------------------------------------
    [zHatTest,~,~] = forward_propagate(x0_TEST,W1,W2,WEnd,b1,b2,bEnd);
    %
    J_Test = sum( sum( 0.5*( z0_TEST - zHatTest ).^2 ) );  
    %
    costVecTest(epochIter) = J_Test / num_ZTest; % Store(TESTING) Cost For SINGLE EPOCH


    %-------------------------------------------------
    % Compute Differences in Errors / Cost Functions
    %            and have it print the info to screen
    %-------------------------------------------------
    if mod( epochIter , hyper_param_vec(6) )==0
        
        %-------------------------------------------------------
        % PRINT INFO TO SCREEN
        %-------------------------------------------------------
        fprintf('#--------------------------------------\n');
        fprintf('    *** Epoch: %d *** \n',epochIter);
        fprintf('Cost (BackP) = %.8f\n',costVec(epochIter) );
        fprintf('Cost (Train) = %.8f\n', J_Train );
        fprintf('Cost  (Test) = %.8f\n\n',J_Test  );
        %pause(7);

    end
    
    
end

fprintf('\n#--------------------------------------\n');
fprintf('   >>>>>> MODEL TRAINING DONE <<<<<< \n');
fprintf('#--------------------------------------');

%-----------------------------
% Data to Give Back to User
%-----------------------------
W1_Save = W1;
W2_Save = W2;
WEnd_Save = WEnd;
%
b1_Save = b1;
b2_Save = b2;
bEnd_Save = bEnd;
%
min_Cost = cost;

