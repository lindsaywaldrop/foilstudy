
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Compute Sobol Indices from aPC Expansion Coefficients!
%
% INPUTS: (1) Coeffs - aPC expansion coefficients vector
%         (2) PolynomialDegree - alpha indices for expansion
%               - each row different expansion coefficient combo for basis
%                 poly's, e.g., row=[1,0,5] --> L_1(x1)*L_0(x2)*L_5(x3)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sFirst,sTot,s2nd,s123] = compute_Sobol_Indices(Coeffs,alphaMAT)

    %-----------------------------------------------
    % N = # of uncertain parameters
    % P = total number of terms in aPC expansion
    %-----------------------------------------------
    [P,N] = size(alphaMAT);

    %--------------------------------------------------------------
    % VAR:  SUM_{j=2}^N Coeffs_j^2 * E[ PSI_j(x1,x2,..,xN)^2] 
    %               (j index in MATLAB notation)
    %--------------------------------------------------------------
    E_PSI_SQR_VEC = compute_Expectation_PSI_Squared(alphaMAT);
    sigmaSqr =  sum( ( Coeffs(2:end).^2 ) .* E_PSI_SQR_VEC(2:end) );


    %-----------------------------
    % Initialization
    %-----------------------------
    sFirst = zeros(1,N);

    %-----------------------------
    % Compute First-Order Index
    %-----------------------------

    for i=1:N % Loop over every uncertain parameter

        %
        % Loop over every row in Coeffs / alpha index matrix for a particular
        %       uncertain parameter, x_i
        ct = 0;
        kVec = [];
        sAux = 0;
        for j=1:P

            %
            % determine if all other indices are zero (or not) in alpha index matrix;
            %        if other indices are not zero, do not use for other first-order index calculation
            use = 1;
            for k=1:N
                if k~=i 
                    if alphaMAT(j,k) ~= 0
                        use = 0;
                    end
                end
            end

            % if only non-zero index corresponds to parameter of interest
            if ( ( alphaMAT(j,i) ~= 0 ) && ( use == 1 ) )

                %poly_index = alphaMAT(j,i);
                %inner_prod = 2 / (2*poly_index + 1); % "2/(2n+1) * delta_mn" if n=m for <Ln,Lm>
                sAux = sAux + Coeffs(j)^2 * E_PSI_SQR_VEC(j);
                kVec = [kVec j];
                ct = ct+1;
            end
        end

        sAux = sAux / sigmaSqr;

        sFirst(i) = sAux;

        %
        % TESTING TO MAKE SURE GRABBING CORRECT ROWS
        %
        % ct
        % kVec
        % for ii=1:ct
        %    alphaMAT(kVec(ii),:) 
        % end

    end


    %-----------------------------
    % Compute TOTAL-Order Index
    %-----------------------------

    for i=1:N % Loop over every uncertain parameter

        %
        % Loop over every row in Coeffs / alpha index matrix for a particular
        %       uncertain parameter, x_i
        ct = 0;
        kVec = [];
        sAux = 0;
        for j=1:P

            % if parameter has non-zero index in alphaMAT
            if ( alphaMAT(j,i) ~= 0 ) 
                sAux = sAux + Coeffs(j)^2 * E_PSI_SQR_VEC(j);
                kVec = [kVec j];
                ct = ct+1;
            end

        end

        sAux = sAux / sigmaSqr;

        sTot(i) = sAux;

        %
        % TESTING TO MAKE SURE GRABBING CORRECT ROWS
        %
        % ct
        % kVec
        % for ii=1:ct
        %    alphaMAT(kVec(ii),:) 
        % end
        % pause();

    end


    %-----------------------------
    % Compute s123... Index
    %-----------------------------

    for i=1:N % Loop over every uncertain parameter

        %
        % Loop over every row in Coeffs / alpha index matrix for a particular
        %       uncertain parameter, x_i
        ct = 0;
        kVec = [];
        sAux = 0;
        for j=1:P

            %
            % determine if all other indices are zero (or not) in alpha index matrix;
            %        if other indices are not zero, do not use for other first-order index calculation
            use = 1;
            for k=1:N
                if alphaMAT(j,k) == 0
                    use = 0;
                end
            end

            % if only non-zero index corresponds to parameter of interest
            if use == 1
                sAux = sAux + Coeffs(j)^2 * E_PSI_SQR_VEC(j);
                kVec = [kVec j];
                ct = ct+1;
            end

        end

        sAux = sAux / sigmaSqr;

        s123(i) = sAux;

        %
        % TESTING TO MAKE SURE GRABBING CORRECT ROWS
        %
        % ct
        % kVec
        % for ii=1:ct
        %    alphaMAT(kVec(ii),:) 
        % end

    end


    %-----------------------------
    % Compute Second-Order Index
    %-----------------------------

    indsVec = 1:N; % indice vector

    for i1=1:N % Loop over every uncertain parameter

        % first index
        ind1 = i1;

        for i2=i1+1:N

            % second index
            ind2 = i2;

            % reinitiate-Auxillary indices vector
            indsAux = zeros( size(indsVec) );

            % store indices being compared and find indices that must be zero
            indsAux(ind1) = ind1;
            indsAux(ind2) = ind2;
            indsTest = find(indsVec-indsAux); % indices we need to be zero

            %
            % Loop over every row in Coeffs / alpha index matrix for a particular
            %       uncertain parameter, x_i
            ct = 0;
            kVec = [];
            sAux = 0;
            for j=1:P

                %
                % determine if all other indices are zero (or not) in alpha index matrix;
                %        if other indices are not zero, do not use for other first-order index calculation
                use = 1;
                for k=1:length(indsTest)
                    if alphaMAT( j,indsTest(k) ) ~= 0
                        use = 0;
                    end
                end

                % if only non-zero indices correspond to parameters of interest
                if ( ( alphaMAT(j,i1) ~= 0 ) && ( alphaMAT(j,i2) ~= 0 ) &&  ( use == 1 ) )
                    sAux = sAux + Coeffs(j)^2 * E_PSI_SQR_VEC(j);
                    kVec = [kVec j];
                    ct = ct+1;
                end



            end

            sAux = sAux / sigmaSqr;

            s2nd(i1,i2) = sAux;

        end

        %
        % TESTING TO MAKE SURE GRABBING CORRECT ROWS
        %
        % ct
        % kVec
        % for ii=1:ct
        %   alphaMAT(kVec(ii),:) 
        % end

    end

