#!/usr/local/Julia/latest/bin/julia
#
using Base.Threads

function generator(SCALE,EdgesPerVertex)
 #   tic()
  N = 2.^SCALE                       # Set  power of number of vertices..

  M = round(Int, EdgesPerVertex .* N)     # Compute total number of edges to generate.

  A = 0.57; B = 0.19;  C = 0.19;   D = 1-(A+B+C)  # Set R-MAT (2x2 Kronecker) coefficients.

  ij = ones(Int, 2, M)         # Initialize index arrays.
  ab = A + B                 # Normalize coefficients.
  c_norm = C/(1 - (A + B))
  a_norm = A/(A + B)
#    begin_time=toq()

 #   tic()
  @threads    for j = 1:M
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
#    loop_time=toq()

    #println("loop Time: " * string(K0time))

 #   tic()
  StartVertex = view(ij,1,:)     # Copy to output. (row vector)
  EndVertex =   view(ij,2,:)       # Copy to output. (row vector)
  #  copy_time=toq()

  #  println("copy output Time: " * string(K0time))

    return StartVertex,EndVertex #, begin_time, loop_time, copy_time
end


SCALE=18
EdgesPerVertex = 16
println("Threads: ",nthreads())
#ut, vt, begin_time, loop_time, copy_time= generator(SCALE, EdgesPerVertex)
ut, vt= generator(SCALE, EdgesPerVertex)
tic()
ut, vt= generator(SCALE, EdgesPerVertex)
K0time=toq()


#println("Begin Time: ", begin_time, " Loop Time: ", loop_time, " Copy time: ", copy_time)
println("Total Time: " * string(K0time))
