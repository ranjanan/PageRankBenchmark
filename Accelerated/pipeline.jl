include("kronecker.jl")
include("write.jl")
function pipeline(scale, edge_factor, num_files)

    Nmax = 2 ^ scale                         
    M = edge_factor * Nmax                
    files = collect(1 : num_files)        

    println("Number of Edges: $M, Maximum Possible Vertex: $Nmax")


    ########################################################
    # Kernel 0: Generate a Graph500 Kronecker graph and save to data files.
    ########################################################

    println("Kernel 0: Generate Graph, Write Edges")

    tic()
    for i in files
        fname1 = "data/K0/from$i.tsv"
        fname2 = "data/K0/to$i.tsv"
        println("  Writing: " * fname1) 
        println("  Writing: " * fname2) 
        srand(i)                                                
        ut, vt = kronecker(scale, edge_factor./ num_files)  
        writeuv(fname1, fname2, ut, vt)
    end
    K0time = toq()

    println("K0 Time: $K0time, Edges/sec: $(M./K0time)")

end
