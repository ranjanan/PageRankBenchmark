module Pagerank

include("kronGraph500NoPerm.jl")
include("io.jl")

function run(scl, EdgesPerVertex)
   info("Scale:", scl)
   info("EdgesPerVertex:", EdgesPerVertex)
   info("Number of threads: ", Threads.nthreads())

   info("Generating data")
   @time ij = kronGraph500NoPerm(scl, EdgesPerVertex)

   filename = "1.tsv"

   info("Writing data:")
   @time writetsv(filename, ij)

   info("Read data")
   @time _ij = readtsv(filename)

   @assert all(ij .== _ij)
   ij
end

end

