%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Model Function
%       
%           eps: SMALL error amount to take into account 
%                   in "model evaluation"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fMAT = Evaluate_Function(DATA,eps,flagScaleOutput,minZ,maxZ)

% Assumes 1 output only
fMAT = zeros(length(DATA(:,1)),1);

% Loop over all input parameter combinations in DATA
for i=1:length(DATA(:,1))

    x1 = DATA(i,1);
    x2 = DATA(i,2);
    x3 = DATA(i,3);
    %x4 = DATA(i,4);
    %x5 = DATA(i,5);
    %x6 = DATA(i,6);
    %x7 = DATA(i,7);

    f1 = (1-x2^2)*cos(2*2*pi*x1);
    f2 = 0.5*sin(2*2*pi*x2);%3*x2^2+0.1;%; %sin(3*2*pi*x2);
    f3 = 2*x3+0.2;
    %f4 = 1.25-x4;
    %f5 = 2*x5^3-3*x5^2+x5+1;
    %f6 = cos(0.5*x6);
    %f7 = (2*x4)*(3*x7)+x7;

    fVal = f1+f2+f3+2;%+f3;%+f4+f5+f6+f7+1;

    fMAT(i,1) = fVal*( 1 + eps*( 2*rand()-1 ) );
end

%-------------------------------------------------------------------------
%                      SCALE OUTPUT DATA TO [0,1] 
%          ** only for evaluating function in error analysis **
%-------------------------------------------------------------------------
if flagScaleOutput
    
    m = 1 / (maxZ-minZ);
    b = -m*minZ;
    fMAT = m * fMAT + b;
    
end 