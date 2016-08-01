using CUDAdrv, CUDAnative

module PageRank
module CUDA

include("kronGraph500NoPerm.jl")
include("../common/common.jl")


#
# Global state
#

type BenchmarkState
    dev::CuDevice
    ctx::CuContext
end

function setup()
    dev = CuDevice(0)
    ctx = CuContext(dev)

    return BenchmarkState(dev, ctx)
end

function teardown(state)
    destroy(state.ctx)
end


#
# Kernel 0
#

type Kernel0Args
    filepath::String
    scl::Int
    EdgesPerVertex::Int
end

kernel0_setup(state, filepath, scl, EdgesPerVertex) = return Kernel0Args(filepath, scl, EdgesPerVertex)

function kernel0(state, args)
   ij1, ij2 = kronGraph500NoPerm_shuffle(args.scl, args.EdgesPerVertex)
   PageRank.write_edges(args.filepath, ij1, ij2)
end

kernel0_teardown(state, args) = return nothing

end
end
