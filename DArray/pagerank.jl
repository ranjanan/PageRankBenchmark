module Pagerank
using DistributedArrays
import DistributedArrays: map_localparts

include("kronGraph500NoPerm.jl")
include("io.jl")

function kernel0(path, scl, EdgesPerVertes)
   info("Generating data")
   @time edges = DArray(I->kronGraph500NoPerm(scl, EdgesPerVertes), (nworkers(), ))

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

