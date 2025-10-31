"""
MyCheb: my implementation of Weideman & Reddy's MATLAB Chebyshev differentiation suite.
"""
module MyCheb 
using LinearAlgebra
"""
mycheb(N):  Computes the Chebyshev nodes of order N.

Based on chebdif from Weideman & Reddy. I started there and just edited out
everything past the part where they figure the Chebyshev points....
leaving me with the following one-liner:

     x = sin(pi*[N-1:-2:1-N]'/(2*(N-1))); % Compute Chebyshev points.

These are not the zeros of ChebyshevT, but of its derivatve..... plus the points +-1.
Note that D_x( T_n(x)) = n U_{n-1}(x)  
where U is Chebyshev U.
So this gives zeros of the ChebyshevU polynomials (as well as +-1, as noted above).
 
ptw 05\02/03
"""
#mycheb(N) = sin(pi*[N-1:-2:1-N]'/(2*(N-1)))
mycheb(N) = sin.(π * collect(N-1:-2:1-N) / (2*(N-1)))

##########################################################################################

"""
    Direct translation to Julia from MATLAB of the code from
    The function [x, DM] = chebdif(N, M) computes the differentiation
    matrices D1, D2, ..., DM on Chebyshev nodes.
    
    Input:
    N:        Size of differentiation matrix.        
    M:        Number of derivatives required (integer).
    Note:     0 < M <= N-1.
    
    Output:
    x:        Chebyshev points
    DM:       DM[:,:,ell] contains ell-th derivative matrix, ell=1..M.
    
    The code implements two strategies for enhanced
    accuracy suggested by W. Don and S. Solomonoff in
    SIAM J. Sci. Comp. Vol. 6, pp. 1253--1268 (1994).
"""
function chebdif(N, M; noflip=false)
    
    #I = Matrix{Float64}(LinearAlgebra.I, N, N)    # Identity matrix
    L = Matrix{Bool}(LinearAlgebra.I, N, N)          # Logical identity matrix
    n1 = floor(Int, N/2)                         # Indices used for flipping trick
    n2 = ceil(Int, N/2)
    k = collect(0:N-1)                           # Compute theta vector
    th = k * π / (N-1)
    x = sin.(π * collect(N-1:-2:1-N) / (2*(N-1))) # Compute Chebyshev points
    
    T = repeat(th/2, 1, N)                       # Create matrix of th/2
    DX = 2 * sin.(T' .+ T) .* sin.(T' .- T)      # Trigonometric identity
    DX = [DX[1:n1,:]; -reverse(reverse(DX[1:n2,:], dims=2), dims=1)]  # Flipping trick
    DX[L] .= 1.0                                 # Put 1's on main diagonal of DX
    
    # Create Toeplitz matrix C
    C = zeros(N, N)
    for i = 1:N, j = 1:N
        C[i,j] = (-1)^(k[i] + k[j])
    end
    C[1,:] *= 2                                  # C is the matrix with
    C[N,:] *= 2                                  # entries c(k)/c(j)  
    C[:,1] /= 2
    C[:,N] /= 2
    
    Z = 1 ./ DX                                  # Z contains entries 1/(x(k)-x(j))
    Z[L] .= 0.0                                  # with zeros on the diagonal
    
    D = Matrix{Float64}(LinearAlgebra.I, N, N)   # D contains diff. matrices
    DM = zeros(N, N, M)                          # Initialize output array
    
    for ell = 1:M
        # Off-diagonals
        D = ell * Z .* (C .* repeat(diag(D), 1, N) - D)
        # Correct main diagonal of D  
        D[L] = -sum(D, dims=2)[:]
        # Store current D in DM
        DM[:,:,ell] = D
    end
    if noflip # The original MATLAB code returns x in reverse order. Ugh.
        return x, DM
    else
        return x[end:-1:1], DM[end:-1:1,end:-1:1,:]
    end
end

end