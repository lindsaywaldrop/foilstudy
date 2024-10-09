
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: create INFORMATION MATRIX 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function INFO_MAT = create_Information_Matrix(N,p,cap_P,param_combo,alphaMAT)

    % Loop over all test points
    for j=1:length( param_combo(:,1) )

        paramVec = param_combo(j,:);

        % Loop over MULTIVARIABLE LEGENDRE POLY
        for i=0:cap_P-1

            % get Multivariable Legendre indices (i+1 bc i starts at 0)
            alphaVec = alphaMAT(i+1,:);

            % compute Multivariable Legendre Poly. at specific point, paramVec
            %         (i+1 bc i starts at 0)
            INFO_MAT(j,i+1) = multi_dim_Legendre_Poly(alphaVec,paramVec);

        end
    end