%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: (1) generates a the Generalized Polynomial Chaos (gPC) 
%               expansion for training data corresponding to an airfoil
%               
%           (2) tests the accuracy of the gPC by:
%                    - recovering the training dataset
%                    - comparing against an independent test dataset
%
%    Author: Nick Battista
%    Date Created: 11/11/2021
%    Modified for Airfoil: April 2025
%    Institution: TCNJ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sCoeffs_cLift, alphaMAT] = perform_gPC_Validation()


%------------------------------------------------------------
% LOADS: gPC Info Data Structures:
%           [1] 'N', 'p', 'cap_P', 'N_subset' (scalars)
%           [2] 'param_combo'
%           [3] 'alphaMAT'
%           [4] 'INFO_MAT'
%------------------------------------------------------------
fprintf('\nLoading gPC Data Structures...\n');
%
load('../results/gPC_Info.mat');
%
fprintf('   --> Done!\n\n');    

%
fprintf('gPC Surrogate Info: \n')
fprintf('    --> # of parameters: %d\n',length(param_combo(1,:)) );
fprintf('    --> Highest polynomial degree: %d\n',p)
fprintf('    --> # of terms in expansion: %d\n',cap_P)



%------------------------------------------------------------------
%   LOAD OUTPUT METRICS:  
%   (** EVALUATED MODEL AT PARTICULAR PARAMETER COMBINATIONS **)
%
%   i.e., loads the training dataset: 
%
%              Training_RESULTS: 263 x 4 matrix 
%                   |--each row different sample
%                   |--col 1: originally sample id
%                   |--col 2: lift data
%                   |--col 3: drag data
%                   |--col 4: center of moment (aerodynamic center)
%       
%              Train_PARAMS: 263 x 3 matrix 
%                   |--each row different parameter combination
%                   |--col 1: Re value
%                   |--col 2: Angle of Attack
%                   |--col 3: Camber
%
%     NOTE: since not all simulations fully ran (ie, converged)
%           the originally sample set will be downsized by the 
%           indices in the first column on Training_RESULTS
%
%------------------------------------------------------------------
fprintf('\nLoading Output Metrics...\n');
%
load('../results/gPC_Training_Dataset.mat');
%
fprintf('    --> Done!\n\n');    


%--------------------------------------------------------------
% MODIFY INFORMATION MATRIX & TRAINING INPUT PARAMETERS 
%                          TO ONLY INCLUDE SIMS THAT FULLY RAN
%--------------------------------------------------------------
SIMS_USED = Training_RESULTS(:,1);
INFO_MAT = INFO_MAT(SIMS_USED,:);
%
fprintf('Simulations fully finished/converged: %d\n',length(SIMS_USED));  


%--------------------------------------------------------------
% DECIDE HOW MUCH TRAINING DATA YOU WANT TO USE!
%           --> most here is 218 (many didn't finish/converge)
%           --> d=6, N=3, (d+1)^N root (param) combinations
%--------------------------------------------------------------
numTrainData = 263; 
%
TRAINING_DATA = Training_RESULTS(1:numTrainData,2:end);
INFO_MAT = INFO_MAT(1:numTrainData,:);
Train_PARAMS = Train_PARAMS(1:numTrainData,:);

% Get metric performance values (from training dataset)
cLift = TRAINING_DATA(1:numTrainData,1);
cDrag = TRAINING_DATA(1:numTrainData,2);
cMoment =  TRAINING_DATA(1:numTrainData,3);


%------------------------------------------------------
% # OF UNKNOWN COEFFS / RANK OF INFOMAT^T * INFOMAT
%------------------------------------------------------
fprintf('# OF UNKNOWN SPECTRAL COEFFS: %d\n',cap_P);
fprintf('RANK OF PSI^T*PSI %d\n', rank( INFO_MAT'*INFO_MAT ) );


%------------------------------------------------------------
% SOLVE FOR COEFFICIENTS VIA LEAST-SQUARES WITH SVD
%          SVD:  A'* A = U Sigma V^*
%------------------------------------------------------------
fprintf('\nComputing gPC Coefficients...\n');
%
% PSEUDO-INVERSE
PS_INV = inv(INFO_MAT'*INFO_MAT)*INFO_MAT';

% Get gPC Coefficients for each output metric
sCoeffs_cDrag = PS_INV*cDrag;
sCoeffs_cLift = PS_INV*cLift;
sCoeffs_cMoment  = PS_INV*cMoment;

fprintf('    --> calculated gPC coeffs!\n');    
fprintf('\n------------------------------------------\n')
%

%------------------------------------------------------------------------
%
%                       RECOVER TRAINING DATASET: 
%
%        --> gPC is able to predict the training data values
%   
%        --> NOTE: same function is used for testing against an 
%                  independent test dataset further down in the code
%
%------------------------------------------------------------------------
fprintf('\n    --> Recover the training dataset\n')
figNum=1;
validate_gPC_Expansion(Train_PARAMS,alphaMAT,...
                                sCoeffs_cDrag,sCoeffs_cLift,sCoeffs_cMoment,...
                                cDrag,cLift,cMoment,figNum);

fprintf('    --> Finished recovering the training dataset\n\n')
fprintf('    --> Hit any [KEY] to continue...\n')
%pause();
fprintf('\n------------------------------------------\n')

                            
%------------------------------------------------------------------------
% LOADS: Performance metrics for the independent test dataset
%        {'Test_PARAMS','Test_RESULTS'}                            
%------------------------------------------------------------------------
load('../results/gPC_Test_Datasets.mat')

Test_PARAMS = Test_PARAMS_noLOG;
Test_RESULTS = Test_RESULTS_noLOG;

LiftTest = Test_RESULTS(:,2);
DragTest = Test_RESULTS(:,3);
CenterMomentTest = Test_RESULTS(:,4);
   

%------------------------------------------------------------------------
%
%                   VALIDATION (PREDICT TEST DATASET)
% 
%    --> gPC is able to predict data for an independent test dataset
%
%------------------------------------------------------------------------
fprintf('\n    --> Validate against the TEST dataset\n')
figNum=3;
validate_gPC_Expansion(Test_PARAMS,alphaMAT,...
                                sCoeffs_cDrag,sCoeffs_cLift,sCoeffs_cMoment,...
                                DragTest,LiftTest,CenterMomentTest,figNum)

                            



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Evaluate Multi-dimensional Legendre Polynomial at specific
%           x,y,z
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = evaluate_MultiDim_Legendre_Poly(sCoeffs,alphaMAT,VEC)


sum = 0;
for i=1:length(alphaMAT(:,1))
    
    prod = 1;
    for j=1:length(alphaMAT(1,:))
        prod = prod * Legendre_Poly(alphaMAT(i,j),VEC(j)); 
    end
        
    sum = sum + sCoeffs(i)*prod;
    
end

val = sum;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: compute multidimension Legendre Polynomial and evaluate at 
%               paramVec
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = multi_dim_Legendre_Poly(alphaVec,paramVec)

val = 1;
for i=1:length(alphaVec)
    
    n = alphaVec(i);  % Get degree with specific single dimensional Leg. Poly  
    x = paramVec(i);  % Get value to evaluate single dimensional Leg. Poly at
    
    % Compute product of all Legendre Poly's evaluated at specific x and of
    %         particular degree
    val = val * Legendre_Poly(n,x);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: get Legendre Polynomial Coefficients: 
%
%       L_n(x) = c0 + c1*x + c2*x^2 +...+ c_n*x^n
%
%       input: n <- degree of Legendre poly
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function coeffs = Legendre_Poly_Coeffs(n)

    if n == 0

        coeffs = 1;      % L_0(x) = 1

    elseif n == 1

        coeffs = [1 0];  % L_1(x) = x

    else

        % Create Legendre Polynomials based off of recurrence relation  
        P_nm1 = 1;
        P_n = [1 0];
        for i=1:(n-1)
            P_np1 = ((2*i+1)*[P_n,0] - i*[0,0,P_nm1])/(i+1); % recurrence
            [P_nm1,P_n] = deal(P_n,P_np1); % shift
        end

        coeffs = P_np1;
        
    end
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: compute nth Legendre polynomial evaluated at x
%           (via Rodrigues' formula)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = Legendre_Poly(n,x)

% from https://en.wikipedia.org/wiki/Legendre_polynomials
if n==0
    val = 1;
elseif n==1
    val = x;
elseif n==2 
    val = 0.5*(3*x^2-1);
elseif n==3
    val = 0.5*(5*x^3-3*x);
elseif n==4
    val = 0.125*( 35*x^4 - 30*x^2 + 3 );
elseif n==5
    val = 0.125*( 63*x^5 - 70*x^3 + 15*x );
elseif n==6 
    val = (0.0625)*( 231*x^6 - 315*x^4 + 105*x^2 - 5 );
elseif n==7
    val = (0.0625)*( 429*x^7 - 693*x^5 + 315*x^3 - 35*x );
elseif n==8
    val = (1/128)*( 6435*x^8 - 12012*x^6 + 6930*x^4 - 1260*x^2 + 35 );
elseif n==9
    val = (1/128)*( 12155*x^9 - 25740*x^7 + 18018*x^5 - 4620*x^3 + 315*x );
elseif n==10
    val = (1/256)*( 46189*x^10 - 109395*x^8 + 90090*x^6 - 30030*x^4 + 3465*x^2 - 63 );

elseif n==11
    val = (1/256)*(x)*( -693 + 15015*x^2 - 90090*x^4 + 218790*x^6 - 230945*x^8 + 88179*x^10 );
    
elseif n==12    
    val = (1/1024)*( 231 - 18018*x^2 + 225225*x^4 - 1021020*x^6 + 2078505*x^8 - 1939938*x^10 + 676039*x^12 );
    
elseif n==13
    val = (1/1024)*(x)*( 3003 - 90090*x^2 + 765765*x^4 - 2771340*x^6 + 4849845*x^8 - 4056234*x^10 + 1300075*x^12 );

elseif n==14
    val = (1/2048)*( -429 + 45045*x^2 - 765765*x^4 + 4849845*x^6 - 14549535*x^8 + 22309287*x^10 - 16900975*x^12 + 5014575*x^14 );

else
    
    % MATLAB's built-in function...gives associated Legendre Polynomial...
    v = legendre(n,x);
    val = v(1);
    
    % RODRIQUEZ FORMULA
    %val = 0;
    %for k=0:n
    %    val = val + nchoosek(n,k)*nchoosek(n+k,k)*( (x-1)/2 )^k;
    %end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: transforms values in [A,B] to desired interval [-1,1]
%
%       TRANSFORM: y=mx+b
%           y(A) = m(A) + b = -1
%           y(B) = m(B) + b = 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function VAL_T = transform_Values_To_Minus1_Plus1(VAL,A,B)


% Setup linear system
MAT = [A 1; B 1];
RHS = [-1;1];
coeffs = MAT \ RHS;

m = coeffs(1);
b = coeffs(2);

VAL_T = m*VAL + b;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: --> validate gPC Expansion against *ACTUAL* Simulated Values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate_gPC_Expansion(REAL_param_MAT,alphaMAT,...
                                sCoeffs_cDrag,sCoeffs_cLift,sCoeffs_cMoment,...
                                cDrag,cLift,cMoment,figNum)


%----------------------------------------------------------
% Loop over ALL parameter combinations model evaluations
%----------------------------------------------------------
for j=1:length( REAL_param_MAT(:,1) )
   
    %-----------------------------------
    % Get Parameter Combination
    %-----------------------------------
    ReTest = REAL_param_MAT(j,1);
    AoATest = REAL_param_MAT(j,2);
    CamberTest = REAL_param_MAT(j,3);
    
    %-----------------------------------------------------
    % Re: Get Transformed FIRST Evaluation Point
    %-----------------------------------------------------
    A = 1e6; B = 1e9;
    Val_Orig = ReTest;
    ReVal = transform_Values_To_Minus1_Plus1(Val_Orig,A,B);

    %-----------------------------------------------------
    % AoA: Get Transformed SECOND Evaluation Point
    %-----------------------------------------------------
    A = 0.0; B = 15.0;
    Val_Orig = AoATest;
    AoAVal = transform_Values_To_Minus1_Plus1(Val_Orig,A,B);
    
    %-----------------------------------------------------
    % Camber: Get Transformed Third Evaluation Point
    %-----------------------------------------------------
    A = 0.005; B = 0.200;
    Val_Orig = CamberTest;
    CamberVal = transform_Values_To_Minus1_Plus1(Val_Orig,A,B);

    %-----------------------------------------------------
    % SPECIFIC PARAMETER COMBINATION
    %-----------------------------------------------------
    VEC = [ReVal AoAVal CamberVal];

    %----------------------------------------------------------------
    %       Evaluate Multidim. Legendre Polys for Each Metric
    %----------------------------------------------------------------
    % Coefficient of Drag
    gPC_CD(j,1) = evaluate_MultiDim_Legendre_Poly(sCoeffs_cDrag,alphaMAT,VEC); 
    %
    % Coefficiet of Lift
    gPC_CL(j,1) = evaluate_MultiDim_Legendre_Poly(sCoeffs_cLift,alphaMAT,VEC); 
    %
    % Center of Moment
    gPC_CM(j,1) = evaluate_MultiDim_Legendre_Poly(sCoeffs_cMoment,alphaMAT,VEC); 
        
end

%--------------------------------------------------------
% PLOT ATTRIBUTES
%--------------------------------------------------------
ms = 50;
ms2= 35;
fs = 20;
colorVec = [0.15 0.70 0.85];  % CAROLINA BLUE!


%--------------------------------------------------------------------
%
%     PLOT gPC Predicted Values against TRUE Simulation Values
%
%--------------------------------------------------------------------

f1 = figure(figNum);
subplot(2,2,1)
plot(1:length(gPC_CD),cDrag,'.','MarkerSize',ms,'Color','k'); hold on;
plot(1:length(gPC_CD),gPC_CD,'.','MarkerSize',ms2,'Color',colorVec); hold on;
leg=legend('Simulated Value','gPC Prediction','Location','NorthOutside');
leg.NumColumns=2;
xlabel('Test Sample ID');
ylabel('Avg. Drag Coefficient');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);
grid on;
grid minor;

subplot(2,2,2)
plot(1:length(gPC_CL),cLift,'.','MarkerSize',ms,'Color','k'); hold on;
plot(1:length(gPC_CL),gPC_CL,'.','MarkerSize',ms2,'Color',colorVec); hold on;
leg=legend('Simulated Value','gPC Prediction','Location','NorthOutside');
leg.NumColumns=2;
xlabel('Test Sample ID');
ylabel('Avg. Lift Coefficient');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);
grid on;
grid minor;

subplot(2,2,3)
plot(1:length(gPC_CL),cLift./cDrag,'.','MarkerSize',ms,'Color','k'); hold on;
plot(1:length(gPC_CL),gPC_CL./gPC_CD,'.','MarkerSize',ms2,'Color',colorVec); hold on;
leg=legend('Simulated Value','gPC Prediction','Location','NorthOutside');
leg.NumColumns=2;
xlabel('n^{th} Training Data Simulation');
ylabel('Avg. Lift/Drag Ratio');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);
grid on;
grid minor;

subplot(2,2,4)
plot(1:length(gPC_CM),cMoment,'.','MarkerSize',ms,'Color','k'); hold on;
plot(1:length(gPC_CM),gPC_CM,'.','MarkerSize',ms2,'Color',colorVec); hold on;
leg=legend('Simulated Value','gPC Prediction','Location','NorthOutside');
leg.NumColumns=2;
xlabel('n^{th} Training Data Simulation');
ylabel('Average cMoment');
set(gca,'FontSize',fs);
set(leg,'FontSize',fs);
grid on;
grid minor;


%-------------------------------------
% Position the Figures Appropriately
%-------------------------------------
f1.Position = [10 10 1075 875]; 


%-----------------------------------------------------------------
% Compute absolute differences btwn gPC predictions & actual data
%-----------------------------------------------------------------
Drag_Diff = gPC_CD-cDrag;
%
Lift_Diff = gPC_CL-cLift;
%
LD_Ratio_Diff = (gPC_CL./gPC_CD) - (cLift./cDrag);
%
Mix_Diff  = gPC_CM-cMoment;


%-----------------------------------------
% Calculate Percent Differences
%-----------------------------------------
Drag_Relative_Error = abs(Drag_Diff) ./ abs(cDrag) * 100;
Lift_Relative_Error = abs(Lift_Diff) ./ abs(cLift) * 100;
LDR_Relative_Error = abs(LD_Ratio_Diff) ./ abs(cLift./cDrag) * 100;
Mix_Relative_Error  = abs(Mix_Diff)  ./ abs(cMoment) * 100;


%----------------------------------------------------
% Plot Percent Difference Information to the Screen
%----------------------------------------------------
fprintf('\n      STATS (Percent Differences):\n');
fprintf('Data | avg  | med  | min  |  max  | stdev\n');
fprintf(' CD  | %.2f | %.2f | %.2f | %.2f | %.2f\n',mean( Drag_Relative_Error ),median(Drag_Relative_Error),min( Drag_Relative_Error ),max( Drag_Relative_Error ),std( Drag_Relative_Error ) );
fprintf(' CL  | %.2f | %.2f | %.2f |  %.2f | %.2f\n',mean( Lift_Relative_Error ),median(Lift_Relative_Error),min( Lift_Relative_Error ),max( Lift_Relative_Error ),std( Lift_Relative_Error ) );
fprintf(' L/D | %.2f | %.2f | %.2f | %.2f | %.2f\n',mean(  LDR_Relative_Error ),median( LDR_Relative_Error),min(  LDR_Relative_Error ),max(  LDR_Relative_Error ),std(  LDR_Relative_Error ) );
fprintf(' CM  | %.2f | %.2f | %.2f | %.2f | %.2f\n\n',mean(  Mix_Relative_Error ), median(Mix_Relative_Error), min( Mix_Relative_Error ), max( Mix_Relative_Error ), std( Mix_Relative_Error ) );


%--------------------------------------------------------------------
%
%    PLOT distributions of each metric's percent differences
%
%--------------------------------------------------------------------


%----------------------------------------
% Histograms: Relative Error: DRAG
%----------------------------------------
f2 = figure(figNum+1);
subplot(4,2,1)
h1a=histogram(Drag_Relative_Error,18,'FaceColor',[0 0.447 0.7410]);
xlabel('Percent Diff: Avg. Drag');
ylabel('PDF');
h1a.Normalization = 'probability';
set(gca,'FontSize',fs);
grid on;
grid minor;

subplot(4,2,2)
h1b=histogram(Drag_Relative_Error,18,'FaceColor',[0.635 0.0780 0.1840]);
xlabel('Percent Diff: Avg. Drag');
ylabel('CDF');
h1b.Normalization = 'cdf';
set(gca,'FontSize',fs);
ylim([0 1.01]);
grid on;
grid minor;


%-------------------------------------------------
% Histograms: Relative Error: LIFT
%-------------------------------------------------
subplot(4,2,3)
h2a=histogram(Lift_Relative_Error,15,'FaceColor',[0 0.447 0.7410]);
xlabel('Percent Diff: Avg. Lift');
h2a.Normalization = 'probability';
ylabel('PDF');
set(gca,'FontSize',fs);
grid on;
grid minor;
%
subplot(4,2,4)
h2b=histogram(Lift_Relative_Error,15,'FaceColor',[0.635 0.0780 0.1840]);
xlabel('Percent Diff: Avg. Lift');
ylabel('CDF');
h2b.Normalization = 'cdf';
set(gca,'FontSize',fs);
ylim([0 1.01]);
grid on;
grid minor;


%-------------------------------------------------
% Histograms: Relative Error: LIFT/DRAG RATIO
%-------------------------------------------------
subplot(4,2,5)
h2a=histogram(LDR_Relative_Error,15,'FaceColor',[0 0.447 0.7410]);
xlabel('Percent Diff: Avg. L/D');
h2a.Normalization = 'probability';
ylabel('PDF');
set(gca,'FontSize',fs);
grid on;
grid minor;
%
subplot(4,2,6)
h2b=histogram(LDR_Relative_Error,15,'FaceColor',[0.635 0.0780 0.1840]);
xlabel('Percent Diff: Avg. L/D');
ylabel('CDF');
h2b.Normalization = 'cdf';
set(gca,'FontSize',fs);
ylim([0 1.01]);
grid on;
grid minor;


%---------------------------------------------------------------------
% Histograms: Relative Error: Center of Moment (Aerodynamic Center)
%---------------------------------------------------------------------
subplot(4,2,7)
h2a=histogram(Mix_Relative_Error,15,'FaceColor',[0 0.447 0.7410]);
xlabel('Percent Diff: Avg. cM');
h2a.Normalization = 'probability';
ylabel('PDF');
set(gca,'FontSize',fs);
grid on;
grid minor;
%
subplot(4,2,8)
h2b=histogram(Mix_Relative_Error,15,'FaceColor',[0.635 0.0780 0.1840]);
xlabel('Percent Diff: Avg. cM');
ylabel('CDF');
h2b.Normalization = 'cdf';
set(gca,'FontSize',fs);
ylim([0 1.01]);
grid on;
grid minor;


%-----------------------------------------
% Position the Plots
%-----------------------------------------
f2.Position = [1000 1 800 2000];








