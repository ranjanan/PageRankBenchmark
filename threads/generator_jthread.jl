#!/usr/local/Julia/latest/bin/julia
#
using Base.Threads

type generator_task
    SCALE::Int
    EdgesPerVertex::Int
    M::Int
    ij1::Array{Int}
    ij2::Array{Int}
end

mytask=generator_task[]


function generator_inner()
    SCALE=18 #mytask[1].SCALE
    EdgesPerVertex = 4 #mytask[1].EdgesPerVertex

  N = 2.^SCALE                       # Set  power of number of vertices..

  np=nthreads()

  r = Base.MersenneTwister(Threads.threadid())

    M = round(Int, EdgesPerVertex .* N)     # Compute total number of edges to generate.

#  M = mytask[1].M

    ij1=Array(Int, M)#mytask[1].ij1
    ij2=Array(Int, M)#mytask[1].ij2

  local_M = div(M,np)
    range = 1:local_M
    lrange=local_M*(threadid()-1) + range

  A = 0.57; B = 0.19;  C = 0.19;   D = 1-(A+B+C)  # Set R-MAT (2x2 Kronecker) coefficients.
  ab = A + B                 # Normalize coefficients.
  c_norm = C/(1 - (A + B))
  a_norm = A/(A + B)

    for j = 1:local_M
       for ib = 1:SCALE            # Loop over each scale.
          k = 1 << (ib-1)
          if rand(r) > ab
              ij1[j] += k
              if rand(r) > c_norm
                  ij2[j] += k
              end
          elseif rand(r) > a_norm
              ij2[j] += k
          end
      end
  end

    return
end

function generator(SCALE,EdgesPerVertex)
    N = 2.^SCALE                       # Set  power of number of vertices..
    M = round(Int, EdgesPerVertex .* N)     # Compute total number of edges to generate.
println(SCALE," ", EdgesPerVertex," " , M)
    push!(mytask, generator_task(SCALE, EdgesPerVertex, M, Array(Int, M), Array(Int, M)))
    ccall(:jl_threading_run, Void, (Any,), Core.svec(generator_inner))
#    return mytask[1].ij1, mytask[1].ij2
end

println("Threads: ",nthreads())

SCALE=18
EdgesPerVertex = 16
tic()
generator(SCALE, EdgesPerVertex)
K0time=toq()

println("Total Time: " * string(K0time))
