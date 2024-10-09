
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: test gPC Expansion for test problem
%
%            -> Qualitative comparison across specific 2D subspace
%            -> Quanitative comparison across that 2D subspace
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_gPC_Expansion(sCoeffs,alphaMAT)

    %------------------------------------------------------------------
    % Compute REAL and gPC predicted value on 2D subspace
    %         (NOTE: Hardcoded function in USER_MODEL_EVALUATION.m)
    %------------------------------------------------------------------
    xVec = linspace(-0.95,0.95,50);
    yVec = linspace(-0.95,0.95,50);
    [X,Y] = meshgrid(xVec,yVec);

    for i=1:length(xVec)
        for j=1:length(yVec)

           % Get specific (x,y,z) value to test in expansion 
           x = xVec(i);
           y = yVec(j);
           z = 0.9;     % Chooses specific 2D subspace

           % Store in vector to pass into expansion
           VEC = [x y z];

           % Evaluate Multidimensional Legendre Polynomial
           Z_Test(j,i) = evaluate_MultiDim_Legendre_Poly(sCoeffs,alphaMAT,VEC); 

           % TEST PROBLEM
           Z_Real(j,i) =  USER_SPECIFIED_MODEL(VEC);

        end 
    end

    %------------------------------------------------
    % PLOT REAL RESPONSE SURFACE VS gPC RESPONSE
    %------------------------------------------------
    f1=figure(1);
    subplot(1,2,1)
    surf(X,Y,Z_Real)
    set(gca,'FontSize',18);
    xlabel('x');
    ylabel('y');
    zlabel('REAL: Response Surface');
    zMin = min(min(Z_Real));
    zMax = max(max(Z_Real));
    h=colorbar;
    caxis([zMin zMax])
    h.Ticks = linspace( zMin, zMax, 5 );
    %
    subplot(1,2,2)
    surf(X,Y,Z_Test)
    set(gca,'FontSize',18);
    xlabel('x');
    ylabel('y');
    zlabel('gPC: Response Surface');
    h=colorbar;
    caxis([zMin zMax])
    h.Ticks = linspace( zMin, zMax, 5 );

    f1.Position = [50 50 1200 400];

    %------------------------------------------------
    % COMPUTE ERROR:
    %    --> Resize data to column vector
    %    --> Compute percent error
    %------------------------------------------------
    len = length(xVec)*length(yVec);
    Z_Real = reshape(Z_Real,len,1);
    Z_Test = reshape(Z_Test,len,1);
    errRel = abs(Z_Real-Z_Test) ./ abs(Z_Real) * 100;
    
    
    fprintf('------------------------------------------\n\n');
    fprintf('gPC Percent Error (on 2D subspace):\n');
    fprintf('   --> AVG:    %.3f\n',mean( errRel)  );
    fprintf('   --> MEDIAN: %.3f\n',median( errRel)  );
    fprintf('   --> MIN:    %.3f\n',min( errRel)  );
    fprintf('   --> MAX:    %.3f\n',max( errRel)  );
    fprintf('   --> STDEV:  %.3f\n\n',std( errRel)  );
    