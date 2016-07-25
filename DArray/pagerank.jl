module Pagerank

include("kronGraph500NoPerm.jl")
include("io.jl")

function run(scl, EdgesPerVertex)
   info("Scale:", scl)
   info("EdgesPerVertex:", EdgesPerVertex)
   info("Number of threads: ", Threads.nthreads())

   info("Generating data")
   @time ij1, ij2 = kronGraph500NoPerm5(scl, EdgesPerVertex)

   filename = "1.tsv"

   info("Writing data:")
   file = BufferedOutputStream(open(filename, "w"))
   @time writetsv(file, ij1, ij2)
   close(file)

   info("Read data")
   file = IOBuffer(Mmap.mmap(open(filename), Vector{UInt8}, (filesize(filename),)))
   @time _ij1, _ij2 = readtsv(file)
   close(file)

   @assert all(ij1 .== _ij1)
   @assert all(ij2 .== _ij2)
   ij1, ij2
end

end

