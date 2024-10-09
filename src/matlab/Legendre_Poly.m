
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: compute nth Legendre polynomial evaluated at x
%
%       --> These polynomials are to be evaluated on [-1,1]
%
%       --> Hardcoded Legendre Polynomials are faster to evaluate than 
%           either Rodriques Formula or MATLAB's Legendre poly function
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

        %----------------------------------------------------------------------------
        % MATLAB's built-in function...gives associated Legendre Polynomial...
        %     --> m = 0 associated Legendre poly is standard Legendre poly.
        %     (*** works but is slower than hardcoding Legendre Polys *** )
        %----------------------------------------------------------------------------
        v = legendre(n,x);
        val = v(1);

        %----------------------------------------------------------------------------
        % RODRIQUES FORMULA 
        %     (*** works but is slower than hardcoding Legendre Polys *** )
        %----------------------------------------------------------------------------
        %val = 0;
        %for k=0:n
        %    val = val + nchoosek(n,k)*nchoosek(n+k,k)*( (x-1)/2 )^k;
        %end

    end