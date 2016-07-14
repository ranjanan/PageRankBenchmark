# execute by running `➜ mpiexec -n $NUM_PROCS julia mpi_pagerank.jl`

import MPI

MPI.Init()

#############
# Utilities #
#############

max_vertex(scale) = 2^scale

total_edges(scale, edges_per_vertex) = edges_per_vertex * max_vertex(scale)

function process_file_range(nfiles, numprocs, rank)
    start_file = ceil(Int, 1 + (nfiles / numprocs) * rank)
    end_file = ceil(Int, (nfiles / numprocs) * (rank + 1))
    return start_file:end_file
end

################################
# Kron Graph Vertex Generation #
################################

function kron_graph_work_arrays(scale, edges_per_vertex)
    m = total_edges(scale, edges_per_vertex)
    start_vertex, end_vertex = Vector{Int}(m), Vector{Int}(m)
    rand_array = Vector{Float64}(2 * scale)
    return start_vertex, end_vertex, rand_array
end

# equivalent to `KronGraph500NoPerm` in the original version
function generate_vertices(scale, edges_per_vertex)
    generate_vertices!(kron_graph_work_arrays(scale, edges_per_vertex)..., scale)
end

function generate_vertices!(start_vertex, end_vertex, rand_array, scale)
    edges = length(start_vertex)

    # set R-MAT (2x2 Kronecker) coefficients
    a = 0.57
    b = 0.19
    c = 0.19
    d = 1.0 - (a + b + c)

    # calculate normalization coefficients
    a_plus_b = a + b
    c_norm = c / (1.0 - a_plus_b)
    a_norm = a / a_plus_b

    # loop over each scale
    @inbounds for j in 1:edges
        rand!(rand_array)
        start_node = one(eltype(start_vertex))
        end_node = one(eltype(end_vertex))
        @inbounds for i in 1:scale
            k = 1 << (i - 2)
            start_bit = rand_array[i] > a_plus_b
            end_bit = rand_array[i + scale] > ifelse(start_bit, c_norm, a_norm)
            start_node += k * start_bit
            end_node += k * end_bit
        end
        start_vertex[j] = start_node
        end_vertex[j] = end_node
    end

    return start_vertex, end_vertex
end


############
# Kernel 0 #
############

# Write vertices to a file. This follows the same TSV specification used in the original
# MATLAB kernel 0 implementation.
function write_vertices!(buffer, file_name, start_vertex, end_vertex)
    for j in eachindex(start_vertex)
        write(buffer, string(start_vertex[j]))
        write(buffer, '\t')
        write(buffer, string(end_vertex[j]))
        write(buffer, '\n')
    end
    open(file_name, "w") do file
        write(file, takebuf_string(buffer))
    end
end

function kernel0(path, files, scale, edges_per_vertex)
    # initialize buffer for edge writing
    buffer = IOBuffer()

    # initialize work arrays
    start_vertex, end_vertex, rand_array = kron_graph_work_arrays(scale, edges_per_vertex)

    # write to files
    for file_index in files
        srand(file_index) # uniquely seed RNG specifically for this file
        file_name = joinpath(path, "$(file_index).tsv")

        # fill vertex vectors with computed edges
        generate_vertices!(start_vertex, end_vertex, rand_array, scale)

        # write vertices to the file
        write_vertices!(buffer, file_name, start_vertex, end_vertex)
    end
end

############
# Kernel 1 #
############

function read_vertices!(start_vertices, end_vertices, file_name, offset)
    new_offset = open(file_name, "r") do file
        i = offset
        while !(eof(file))
            i += 1
            start_vertices[i] = parse(Int, readuntil(file, '\t'))
            end_vertices[i] = parse(Int, readuntil(file, '\n'))
        end
        return i
    end
    return new_offset
end

function kernel1_read(path, files, scale, edges_per_vertex)
    nvertices = total_edges(scale, edges_per_vertex) * length(files)
    start_vertices = Vector{Int}(nvertices)
    end_vertices = Vector{Int}(nvertices)
    offset = 0
    for file_index in files
        file_name = joinpath(path, "$(file_index).tsv")
        offset = read_vertices!(start_vertices, end_vertices, file_name, offset)
    end
    perm = sortperm(start_vertices)
    return start_vertices[perm], end_vertices[perm]
end

edge_bound_indices(numprocs, nedges) = div((0:(numprocs - 1)) * nedges, numprocs) + 1

function kernel1_edge_bounds(numprocs, start_vertices, scale)
    # extract edge boundaries from this process's starting vertices
    start_bounds = start_vertices[edge_bound_indices(numprocs, length(start_vertices))]

    # ensure that there's no overlaps in our selected edge boundaries
    start_bounds[1] = 1
    for i in 2:length(start_bounds)
        current, previous = start_bounds[i], start_bounds[i - 1]
        if current <= previous
            start_bounds[i] = previous + 1
        end
    end

    # from start_bounds, compute the edge boundaries for the ending vertices
    end_bounds = vcat(start_bounds[2:end]-1, max_vertex(scale))

    return start_bounds, end_bounds
end

# function kernel1_sample_sort(start_vertices, end_vertices, scale)
#     numprocs = MPI.Comm_size(MPI.COMM_WORLD)
#     myrank = MPI.Comm_rank(MPI.COMM_WORLD)
#
#     start_bounds, end_bounds = kernel1_edge_bounds(numprocs, start_vertices, scale)
#
#     # prepare reception of the correct edges from other processes
#     for rank in 0:(numprocs - 1)
#         i = rank + 1
#         lower = lower_bounds[i]
#         upper = upper_bounds[i]
#         inbounds = find(x -> lower <= x <= upper, start_vertices)
#         # allocate reception buffer for start_vertices
#         MPI.Irecv!(start_vertices[inbounds], rank, rank, MPI.COMM_WORLD)
#         # allocate reception buffer for end_vertices
#         MPI.Irecv!(end_vertices[inbounds], rank, i, MPI.COMM_WORLD)
#     end
#
#     # send this process's edges to the correct processes
#     for rank in 0:(numprocs - 1)
#         i = rank + 1
#         lower = lower_bounds[i]
#         upper = upper_bounds[i]
#         inbounds = find(x -> lower <= x <= upper, start_vertices)
#         # send start_vertices
#         MPI.Isend(start_vertices[inbounds], rank, rank, MPI.COMM_WORLD)
#         # send end_vertices
#         MPI.Isend(end_vertices[inbounds], rank, i, MPI.COMM_WORLD)
#     end
# end

##########
# main() #
##########

function main()
    # process info
    is_master_process = MPI.Comm_rank(MPI.COMM_WORLD) == 0
    numprocs = MPI.Comm_size(MPI.COMM_WORLD)
    myrank = MPI.Comm_rank(MPI.COMM_WORLD)

    # problem parameters
    scale = 18
    edges_per_vertex = 16     # must be power of 2
    nfiles = edges_per_vertex # must be power of 2
    files = process_file_range(nfiles, numprocs, myrank)

    # set up directoriess
    data_path = joinpath(pwd(), "pagerankdata")
    kernel0_path = joinpath(data_path, "kernel0")

    if is_master_process
        !(isdir(data_path)) && mkdir(data_path)
        !(isdir(kernel0_path)) && mkdir(kernel0_path)
    end

    # kernel 0
    is_master_process && print("running kernel 0...")
    is_master_process && tic()

    MPI.Barrier(MPI.COMM_WORLD)
    kernel0(kernel0_path, files, scale, edges_per_vertex)
    MPI.Barrier(MPI.COMM_WORLD)

    is_master_process && println("done (took $(toq()) seconds)")

    # kernel 1
    is_master_process && print("running kernel 1...")
    is_master_process && tic()

    MPI.Barrier(MPI.COMM_WORLD)
    # kernel1()
    MPI.Barrier(MPI.COMM_WORLD)

    is_master_process && println("done (took $(toq()) seconds)")
end

main()

MPI.Finalize()
