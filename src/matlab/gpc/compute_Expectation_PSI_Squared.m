
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: compute the Expectation of PSI^2 and store each into vector
%                          ( E[ PSI_j^2(x) ] )
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E_PSI_SQR_VEC = compute_Expectation_PSI_Squared(alphaMAT)

    % Loop over each term in the gPC expansion
    for i=1:length(alphaMAT(:,1))

        % Loop over each index in alphaMAT 
        %       (loops over # of uncertain parameters)
        prod = 1;
        for j=1:length(alphaMAT(1,:))

            prod = prod * 1 / (2*alphaMAT(i,j)+1);

        end

        E_PSI_SQR_VEC(i,1) = prod;

    end