module Pagerank
using DistributedArrays
import DistributedArrays: map_localparts

include("kronGraph500NoPerm.jl")
include("io.jl")

function kernel0(path, scl, EdgesPerVertes)
   n = 2^scl
   m = EdgesPerVertex * n # Total number of vertices

   # Make sure that we distribute the workload over the workers.
   EdgesPerWorker = m ÷ nworkers()
   surplus = m % nworkers()

   info("Generating data")
   lastWorker = maximum(workers())
   @time begin
      rrefs = map(workers()) do id
         nEdges = EdgesPerWorker
         nEdges += ifelse(id == lastWorker, surplus, 0)
         kronGraph500NoPerm(scl, nEdges)
      end
      edges = DArray(rrefs)
   end

   filenames = Dict(id => joinpath(path, "$i.tsv") for (id, i) in zip(workers(), 1:nworkers()))
   info("Writing data:")
   @time map_localparts(edges) do local_edges
      writecsv(filename, filenames[myid()], local_edges)
   end
   return values(filenames)
end

function kernel1(filenames)
   info("Read data")
   @time begin
      rrefs = map(zip(filennames, workers())) do iter
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

   info("Executing kernel 0")
   @time filenames = kernel0(path, scl, EdgesPerVertex)

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

