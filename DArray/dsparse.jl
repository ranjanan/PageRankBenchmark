# Input is a vector of futures
function create_adj_matrix(fvectors)      
   rrefs = map(fvectors) do fvec
      remotecall(fvec.where) do
         edges = fetch(fvec)
         # SparseMatrixCSC only accepts I, J, V for construction
         # So this is fairly expensive
         I = Vector{Int64}(length(edges))
         J = Vector{Int64}(length(edges))
         V = Vector{Int64}(length(edges))

         min_i = first(edges)[1] - 1
         max_j = 0
         max_i = 0

         for ind in eachindex(edges, I, J, V)
            i, j = edges[ind]
            i = i - min_i # localparts of sparse matrix need to start at 1
            max_j = ifelse(j > max_j, j, max_j)
            max_i = ifelse(i > max_i, i, max_i)
            I[ind] = i
            J[ind] = j
            V[ind] = 1
         end
         (I, J, V, max_i, max_j)
      end
   end

   # Collect the maxium
   max_js = map(rrefs) do rref
      remotecall_fetch(r -> fetch(r)[5], rref.where, rref)
   end

   max_j = maximum(max_js)

   # Construct sparse array
   lparts = map(rrefs) do rref
      remotecall(rref.where, rref) do ref
         (I, J, V, max_i, mj) = fetch(ref) # max_i needs to be local, max_j needs to be global
         sparse(I, J, V, max_i, max_j)
      end
   end
   return DArray(reshape(lparts, (length(lparts), 1)))
end
