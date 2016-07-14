#!/usr/local/Julia/latest/bin/julia
#
using Base.Threads


function generator_inner()
SCALE=18
EdgesPerVertex = 4 #16

   # tic()

  N = 2.^SCALE                       # Set  power of number of vertices..

  np=nthreads()
  M = round(Int, EdgesPerVertex .* N)     # Compute total number of edges to generate.

  M = div(M,np)

  A = 0.57; B = 0.19;  C = 0.19;   D = 1-(A+B+C)  # Set R-MAT (2x2 Kronecker) coefficients.
  ij = ones(Int, 2, M)         # Initialize index arrays.

  ab = A + B                 # Normalize coefficients.
  c_norm = C/(1 - (A + B))
  a_norm = A/(A + B)

   # begin_time=toq()

  #  tic()


    for j = 1:M
      for ib = 1:SCALE            # Loop over each scale.
          k = 1 << (ib-1)
          if rand() > ab
              ij[1,j] += k
              if rand() > c_norm
                  ij[2,j] += k
              end
          elseif rand() > a_norm
              ij[2,j] += k
          end
      end
  end

 #   loop_time=toq()
 #   tic()
  StartVertex = view(ij,1,:)     # Copy to output. (row vector)
  EndVertex =   view(ij,2,:)       # Copy to output. (row vector)
#    copy_time=toq()
    return

end

println("Threads: ",nthreads())
#generator_inner()
tic()
#generator_inner()
ccall(:jl_threading_run, Void, (Any,), Core.svec(generator_inner))
K0time=toq()


#println("Begin Time: ", begin_time, " Loop Time: ", loop_time, " Copy time: ", copy_time)
println("Total Time: " * string(K0time))
