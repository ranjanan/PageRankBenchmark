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

   info("Writing data:")
   @time map_localparts(edges) do local_edges
      filename = joinpath(path, "$(myid()).tsv")
      writecsv(filename, local_edges)
   end
end

function kernel1(path)
   info("Read data")
   rrefs = map(workers()) do id
      filename = joinpath(path, "$id.tsv")
      remotecall(readtsv, id, filename)
   end
   edges = DArray(rrefs)
end

function run(path, scl, EdgesPerVertex)
   info("Scale:", scl)
   info("EdgesPerVertex:", EdgesPerVertex)
   info("Number of workers: ", nworkers())

   info("Executing kernel 0")
   @time kernel0(path, scl, EdgesPerVertex)

   info("Executing kernel 1")
   @time kernel1(path)
end

end

