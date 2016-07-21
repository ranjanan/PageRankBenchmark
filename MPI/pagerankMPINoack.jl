# using MPI

# MPI.Init()

# if id == 0
#     a = [22,7,13,18,2,17,1,14]
# elseif id == 1
#     a = [20,6,10,24,15,9,21,3]
# elseif id == 2
#     a = [16,19,23,4,11,12,5,8]
# else
#     error("only three processes supported in this example")
# end

function sortSample(a)

    @show ct = MPI.COMM_WORLD
    @show id = MPI.Comm_rank(ct)
    @show p  = MPI.Comm_size(ct)

    as = sort(a)
    @show as

    n = length(a)

    @show splitters_local = [as[round(Int, (n + 1)/p*i)] for i = 1:(p - 1)]

    @show splitters_all = MPI.Gather(splitters_local, 0, ct)

    @show splitters_final = Array(Int, 2)

    if id == 0
        println("What?")
        sort!(splitters_all)
        println("Noack")
        for i = 1:p - 1
            splitters_final[i] = splitters_all[round(Int, (length(splitters_all) + 1)/p*i)]
        end
        println("Jensen")
    end

    @show splitters_final

    MPI.Bcast!(splitters_final, 0, ct)

    @show offsets_local = Cint[0]
    i = j = 1
    while j <= length(splitters_final)
        if as[i] >= splitters_final[j]
            push!(offsets_local, i - 1)
            j += 1
        end
        i += 1
    end
    push!(offsets_local, n)
    @show scounts = diff(offsets_local)

    @show scounts_global = MPI.Allgather(scounts, ct)

    @show splitters_final
    @show id, scounts
    @show scounts_global

    rcounts = zeros(Cint, length(scounts))
    for i = 0:p - 1
        rcounts[i + 1] += scounts_global[i*length(scounts) + id + 1]
    end
    # @show id, rcounts

    b = sort!(MPI.Alltoallv(as, scounts, rcounts, ct))

    return b
end
    # @show id, 

# MPI.Finalize()
