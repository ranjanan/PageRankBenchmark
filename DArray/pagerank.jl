module Pagerank
using DistributedArrays

include("kronGraph500NoPerm.jl")
include("io.jl")

function kernel0(filenames, scl, EdgesPerVertex)
   n = 2^scl
   m = EdgesPerVertex * n # Total number of vertices

   # Make sure that we distribute the workload over the workers.
   EdgesPerWorker = m ÷ nworkers()
   surplus = m % nworkers()

   @assert length(filenames) == nworkers()

   lastWorker = maximum(workers())
   @sync begin
      for (id, filename) in zip(workers(), filenames)
         nEdges = EdgesPerWorker
         nEdges += ifelse(id == lastWorker, surplus, 0)
         @async remotecall_wait(kronGraph500, id, filename, scl, nEdges)
      end
   end
end

function kernel1(filenames)
   info("Read data")
   @time begin
      rrefs = map(zip(filenames, workers())) do iter
         filename, id = iter
         remotecall(readtsv, id, filename)
      end
      edges = DArray(rrefs)
   end

   info("Sort edges")
   @time sorted_edges = sort(edges)
   close(edges)

   info("Turn into adjacency matrix")
   @time begin
      rrefs = map(workers()) do id
         remotecall(id) do
            ledges = localpart(sorted_edges)
            # SparseMatrixCSC only accepts I, J, V for construction
            # So this is fairly expensive
            I = Vector{Int64}(length(ledges))
            J = Vector{Int64}(length(ledges))
            V = Vector{Int64}(length(ledges))

            min_i = first(ledges)[1] - 1
            max_j = 0
            max_i = 0

            for ind in eachindex(ledges, I, J, V)
               i, j = ledges[ind]
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
      adj_matrix = DArray(reshape(lparts, (length(lparts), 1)))
   end
end

function run(path, scl, EdgesPerVertex)
   info("Scale:", scl)
   info("EdgesPerVertex:", EdgesPerVertex)
   info("Number of workers: ", nworkers())

   filenames = collect(joinpath(path, "$i.tsv") for i in 1:nworkers())

   info("Executing kernel 0")
   @time kernel0(filenames, scl, EdgesPerVertex)

   # Shuffle the filenames so that we minimise cache effect
   # TODO ideally we would like to make sure that no processor reads in
   # its own file.
   shuffle!(filenames)

   info("Executing kernel 1")
   @time kernel1(filenames)
end

# Two helper functions to make sort work on tuples
Base.typemin{T}(::Type{Tuple{T,T}}) = (typemin(T),typemin(T))
Base.typemax{T}(::Type{Tuple{T,T}}) = (typemax(T),typemax(T))

end

