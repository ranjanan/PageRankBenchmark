using CUDAdrv

using BenchmarkTools

dev = CuDevice(0)
ctx = CuContext(dev)

include("kronGraph500NoPerm.jl")
println(@benchmark PageRank.kronGraph500NoPerm(20,16) evals=1 seconds=5)

destroy(ctx)

include("../DArray/kronGraph500NoPerm.jl")
println(@benchmark PageRank.kronGraph500NoPerm4(20,16) evals=1 seconds=5)
