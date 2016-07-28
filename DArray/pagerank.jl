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
   @time sort(edges)
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

