function x = mycheb(N)

%  Computes the Chebyshev nodes of order N.
%
% Based on chebdif from Weideman & Reddy. I started there and just edited out
% everything past the part where they figure the Chebyshev points....
% leaving me with the following one-liner:

     x = sin(pi*[N-1:-2:1-N]'/(2*(N-1))); % Compute Chebyshev points.

% These are not the zeros of ChebyshevT, but of its derivatve..... plus the points +-1.
% Note that D_x( T_n(x)) = n U_{n-1}(x)  
% where U is Chebyshev U.
% So this gives zeros of the ChebyshevU polynomials (as well as +-1, as noted above).
% 
% ptw 05\02/03

