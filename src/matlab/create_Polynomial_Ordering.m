

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: get convention for how to setup multivariable Legendre
%           polynomial
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function alphaMAT = create_Polynomial_Ordering(N,p)

    q = N;  % # of varies parameters (qualities)
    d = p;  % highest degree basis polynomial

    %-----------------------------------------------------------
    % Start with column vector of poly orders {0,1,2,...,d}
    %-----------------------------------------------------------
    alphaMAT = (0:1:d)';

    %----------------------------------------------------------------
    % Loop over how many qualities we're varying 
    %(minus 1, since starting with  a single alpha vector already)
    %----------------------------------------------------------------
    for n=1:q-1

        %-----------------------------------------------------
        % empty vector for 'new' leftmost vector in alphaMAT
        %-----------------------------------------------------
        leftAdd = [];

        %---------------------------------------------------------------
        % CREATE left-most vector --> loop over each plausible degree
        %---------------------------------------------------------------
        for j=0:d

            %------------------------------------------------------------------------
            % create a vector of all 'j' values of length of previous alpha matrix
            %------------------------------------------------------------------------
            leftAux = j*ones( length( alphaMAT(:,1) ) , 1 );

            %-----------------------------------------------------
            % concatenate new left-most vector for alphaMAT
            %-----------------------------------------------------
            leftAdd = [leftAdd; leftAux];

        end

        %--------------------------------------------------------------
        % create empty array to stack alphaMAT in a certain # of times
        %--------------------------------------------------------------
        alphaStack = [];

        %-----------------------------------------------------
        % STACK the alphaMAT we already have!
        %-----------------------------------------------------
        for j=0:d


            %------------------------------------------------
            % stack the pre-existing alphaMAT from before
            %------------------------------------------------
            alphaStack = [alphaStack; alphaMAT];

        end

        %----------------------------------------------------------------------
        % REDEFINE alphaMAT to include new leftVector and stacked alphaMat!
        %----------------------------------------------------------------------
        alphaMAT = [leftAdd alphaStack];

    end

    %------------------------------------------------------
    % CUT OUT COMBINATIONS WHERE SUM INDICES > d
    %------------------------------------------------------
    alphaNEW = [];
    for i=1:length(alphaMAT(:,1))

        if sum( alphaMAT(i,:) , 2) <= p
            alphaNEW = [alphaNEW; alphaMAT(i,:)];
        end

    end

    alphaMAT = alphaNEW;

    %------------------------------------------------------
    % ORIGINAL HARDCODED METHOD
    %------------------------------------------------------
    % ct=0;
    % for i=0:p
    %     for j=0:p
    %         for k=0:p
    %            % if sum of powers is <= p, store index combo! 
    %            if i+j+k<=p
    %                ct = ct+1;
    %                alphaMAT(ct,1) = i;
    %                alphaMAT(ct,2) = j;
    %                alphaMAT(ct,3) = k;
    %            end
    %         end
    %     end
    % end