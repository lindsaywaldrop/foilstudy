
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: subsample all parameter combinations by taking the 
%           N_subset closest number to origin
%           (parameters are all scaled to [-1,1] here)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function param_combo = sample_Parameter_Combinations(N_subset,param_combo)


    for i=1:length(param_combo(:,1))

       % compute overall distance for origin, with points in particular row
       %        of parameter combination
       sum = 0;
       for j=1:length( param_combo(1,:))
           sum = sum + param_combo(i,j)^2;
       end

       % distance from origin
       distanceVec(i) = sum;

    end

    %----------------------------------------------------------------------
    % Re-order the data in distanceVec, get orig inds, take first N_subset
    %----------------------------------------------------------------------
    [~,inds] = sort( distanceVec );
    %
    param_combo = param_combo(inds(1:N_subset),:);
