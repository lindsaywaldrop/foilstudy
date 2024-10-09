
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: combine all plausible collocation points, i.e.,
%           combine roots of Legendre polynomial --> uses same machinery
%           as creating the alphaMAT indices of the polynomial.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function param_combo = compute_All_Collo_PT_Combinations(N,poly_roots)


    %-----------------------------------------------------------
    % Start with column vector of POLY ROOTS
    %-----------------------------------------------------------
    rootMAT = poly_roots;
    numRoots = length(poly_roots); 

    %----------------------------------------------------------------
    % Loop over how many qualities we're varying 
    %(minus 1, since starting with  a root vector already)
    %----------------------------------------------------------------
    for n=1:N-1

        %-----------------------------------------------------
        % empty vector for 'new' leftmost vector in rootMAT
        %-----------------------------------------------------
        leftAdd = [];

        %---------------------------------------------------------------
        % CREATE left-most vector --> loop over each plausible degree
        %---------------------------------------------------------------
        for j=1:numRoots

            %------------------------------------------------------------------------
            % create a vector of all 'j' values of length of previous root matrix
            %------------------------------------------------------------------------
            leftAux = poly_roots(j)*ones( length( rootMAT(:,1) ) , 1 );

            %-----------------------------------------------------
            % concatenate new left-most vector for rootMAT
            %-----------------------------------------------------
            leftAdd = [leftAdd; leftAux];

        end

        %--------------------------------------------------------------
        % create empty array to stack rootMAT in a certain # of times
        %--------------------------------------------------------------
        rootStack = [];

        %-----------------------------------------------------
        % STACK the rootMAT we already have!
        %-----------------------------------------------------
        for j=1:numRoots


            %------------------------------------------------
            % stack the pre-existing rootMAT from before
            %------------------------------------------------
            rootStack = [rootStack; rootMAT];

        end

        %----------------------------------------------------------------------
        % REDEFINE rootMAT to include new leftVector and stacked rootMAT!
        %----------------------------------------------------------------------
        rootMAT = [leftAdd rootStack];

    end

    %------------------------------------------------
    % DEFINE param_combo for output!
    %------------------------------------------------
    param_combo = rootMAT;


    %------------------------------------------------
    % ORIGINAL METHODOLOGY
    %------------------------------------------------
    % % allocate memory
    % param_combo = zeros(length(poly_roots)^N,N);
    % 
    % % compute ALL parameter combinations
    % ct = 0;
    % for i=1:length(poly_roots)
    %     for j=1:length(poly_roots)
    %         for k=1:length(poly_roots)
    %             ct = ct+1;
    %             param_combo(ct,1) = poly_roots(i);
    %             param_combo(ct,2) = poly_roots(j);
    %             param_combo(ct,3) = poly_roots(k);
    %         end
    %     end
    % end
    % size(param_combo)
