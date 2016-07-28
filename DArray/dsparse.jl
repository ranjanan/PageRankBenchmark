# Input is a vector of futures
function create_adj_matrix(fvectors, maxI, maxJ)
   firstWorker = first(workers())
   lastWorker = last(workers())

   rrefs = map(fvectors) do fvec
      remotecall(fvec.where) do
         edges = fetch(fvec)
         # SparseMatrixCSC only accepts I, J, V for construction
         # So this is fairly expensive
         I = Vector{Int64}(length(edges))
         J = Vector{Int64}(length(edges))
         V = Vector{Int64}(length(edges))

         min_i = ifelse(myid() == firstWorker, 1, first(edges)[1]) - 1
         max_i = last(edges)[1]

         for ind in eachindex(edges, I, J, V)
            i, j = edges[ind]
            i = i - min_i # localparts of sparse matrix need to start at 1
            I[ind] = i
            J[ind] = j
            V[ind] = 1
         end
         (I, J, V, max_i, min_i)
      end
   end

   # Collect the maxium
   max_is = map(rrefs) do rref
      remotecall_fetch(r -> fetch(r)[4], rref.where, rref)
   end

   surplus = maxI - maximum(max_is)

   # Construct sparse array
   lparts = map(rrefs) do rref
      remotecall(rref.where, rref) do ref
         (I, J, V, max_i, min_i) = fetch(ref) # max_i needs to be local, maxJ needs to be global
         max_i = ifelse(myid() == lastWorker, max_i + surplus, max_i) - min_i
         sparse(I, J, V, max_i, maxJ)
      end
   end
   return DArray(reshape(lparts, (length(lparts), 1)))
end
