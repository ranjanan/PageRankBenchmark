using CUDAdrv, CUDAnative

module PageRank

include("kronGraph500NoPerm.jl")
include("../io/io.jl")

function kernel0(filename, scl, EdgesPerVertex)
   n = 2^scl
   m = EdgesPerVertex * n # Total number of vertices

   ij1, ij2 = kronGraph500NoPerm_shuffle(scl, EdgesPerVertex)
   PagerankIO.write_edges(filename, ij1, ij2)
end

function run(path, scl, EdgesPerVertex)
   info("Scale:", scl)
   info("EdgesPerVertex:", EdgesPerVertex)

   info("Executing kernel 0")
   @time kernel0(joinpath(path, "kernel0.tsv"), scl, EdgesPerVertex)
end

end
