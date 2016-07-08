module PageRank

# Almost literal translation of MATLAB code
function kronGraph500NoPerm(scl, EdgesPerVertex)
# Graph500NoPerm: Generates graph edges using the same 2x2 Kronecker algorithm (R-MAT) as the Graph500 benchmark, but no permutation of vertex labels is performed.
# IO user function.
#   Usage:
#     StartVertex EndVertex = Graph500NoPerm(scl, edgefactor)
#   Inputs:
#     scl = integer scale factor that sets the max number of vertices to 2^scl
#     EdgesPerVertex = sets the total number of edges to M = K*N;
#  Outputs:
#     StartVertex = M vector of integer start vertices in the range [1,N]
#     EndVertex = M vector of integer end vertices in the range [1,N]

    n = 2^scl                          # Set  power of number of vertices..

    m = EdgesPerVertex * n             # Compute total number of edges to generate.

    a, b, c = 0.57, 0.19, 0.19
    d = 1 - (a + b + c)                # Set R-MAT (2x2 Kronecker) coefficeints.

    ij = ones(2, m)                    # Initialize index arrays.
    ab = a + b                         # Normalize coefficients.
    c_norm = c/(1 - (a + b))
    a_norm = a/(a + b)

    for ib = 1:scl                     # Loop over each scale.
        ii_bit = rand(m) .> ab
        jj_bit = rand(m) .> ( c_norm * ii_bit + a_norm * !(ii_bit) )
        ij = ij + 2^(ib - 1) * [ii_bit'; jj_bit']
    end

    return ij[1,:], ij[2,:]
end

# Devectorize
function kronGraph500NoPerm2(scl, EdgesPerVertex)
# Graph500NoPerm: Generates graph edges using the same 2x2 Kronecker algorithm (R-MAT) as the Graph500 benchmark, but no permutation of vertex labels is performed.
# IO user function.
#   Usage:
#     StartVertex EndVertex = Graph500NoPerm(scl, edgefactor)
#   Inputs:
#     scl = integer scale factor that sets the max number of vertices to 2^scl
#     EdgesPerVertex = sets the total number of edges to M = K*N;
#  Outputs:
#     StartVertex = M vector of integer start vertices in the range [1,N]
#     EndVertex = M vector of integer end vertices in the range [1,N]

    n = 2^scl                             # Set  power of number of vertices..

    m = EdgesPerVertex * n                # Compute total number of edges to generate.

    a, b, c = 0.57, 0.19, 0.19
    d = 1 - (a + b + c)                   # Set R-MAT (2x2 Kronecker) coefficeints.

    ij1, ij2 = ones(Int, m), ones(Int, m) # Initialize index arrays.
    ab = a + b                            # Normalize coefficients.
    c_norm = c/(1 - (a + b))
    a_norm = a/(a + b)

    for ib = -1:scl - 2                   # Loop over each scale.
        @inbounds for i = 1:m
            ii_bit  = rand() > ab
            jj_bit  = rand() > ifelse(ii_bit, c_norm, a_norm)
            sc      = 1 << ib
            ij1[i] += sc * ii_bit
            ij2[i] += sc * jj_bit
        end
    end

    return ij1, ij2
end

# Move rng out of inner loop
function kronGraph500NoPerm3(scl, EdgesPerVertex)
# Graph500NoPerm: Generates graph edges using the same 2x2 Kronecker algorithm (R-MAT) as the Graph500 benchmark, but no permutation of vertex labels is performed.
# IO user function.
#   Usage:
#     StartVertex EndVertex = Graph500NoPerm(scl, edgefactor)
#   Inputs:
#     scl = integer scale factor that sets the max number of vertices to 2^scl
#     EdgesPerVertex = sets the total number of edges to M = K*N;
#  Outputs:
#     StartVertex = M vector of integer start vertices in the range [1,N]
#     EndVertex = M vector of integer end vertices in the range [1,N]

    n = 2^scl                             # Set  power of number of vertices..

    m = EdgesPerVertex * n                # Compute total number of edges to generate.

    a, b, c = 0.57, 0.19, 0.19
    d = 1 - (a + b + c)                   # Set R-MAT (2x2 Kronecker) coefficeints.

    ij1, ij2 = ones(Int, m), ones(Int, m) # Initialize index arrays.
    ab = a + b                            # Normalize coefficients.
    c_norm = c/(1 - (a + b))
    a_norm = a/(a + b)

    rand1buf = Array(Float64, m)
    rand2buf = Array(Float64, m)
    for ib = -1:scl - 2                   # Loop over each scale.
        rand!(rand1buf)
        rand!(rand2buf)
        sc = 1 << ib
        @inbounds for i = 1:m
            ii_bit  = rand1buf[i] > ab
            jj_bit  = rand2buf[i] > ifelse(ii_bit, c_norm, a_norm)
            ij1[i] += sc * ii_bit
            ij2[i] += sc * jj_bit
        end
    end

    return ij1, ij2
end

# FASTEST VERSION
# change loop order to reduce loads and stores
function kronGraph500NoPerm4(scl, EdgesPerVertex)
# Graph500NoPerm: Generates graph edges using the same 2x2 Kronecker algorithm (R-MAT) as the Graph500 benchmark, but no permutation of vertex labels is performed.
# IO user function.
#   Usage:
#     StartVertex EndVertex = Graph500NoPerm(scl, edgefactor)
#   Inputs:
#     scl = integer scale factor that sets the max number of vertices to 2^scl
#     EdgesPerVertex = sets the total number of edges to M = K*N;
#  Outputs:
#     StartVertex = M vector of integer start vertices in the range [1,N]
#     EndVertex = M vector of integer end vertices in the range [1,N]

    n = 2^scl                             # Set  power of number of vertices..

    m = EdgesPerVertex * n                # Compute total number of edges to generate.

    a, b, c = 0.57, 0.19, 0.19
    d = 1 - (a + b + c)                   # Set R-MAT (2x2 Kronecker) coefficeints.

    ij1, ij2 = ones(Int, m), ones(Int, m) # Initialize index arrays.
    ab = a + b                            # Normalize coefficients.
    c_norm = c/(1 - (a + b))
    a_norm = a/(a + b)

    randbuf = Array(Float64, 2scl)
    @inbounds for i = 1:m
        rand!(randbuf)
        ij1_i = zero(eltype(ij1))
        ij2_i = zero(eltype(ij2))
        @inbounds for ib = 1:scl                   # Loop over each scale.
            sc = 1 << (ib - 2)
            ii_bit  = randbuf[ib] > ab
            jj_bit  = randbuf[ib+scl] > ifelse(ii_bit, c_norm, a_norm)
            ij1_i += sc * ii_bit
            ij2_i += sc * jj_bit
        end
        ij1[i] = ij1_i
        ij2[i] = ij2_i
    end

    return ij1, ij2
end

end