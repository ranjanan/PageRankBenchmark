#
@everywhere include("KronGraph500NoPerm.jl")
@everywhere include("StrFileWrite.jl")
@everywhere include("StrFileRead.jl")

function PageRankPipeline(SCALE,EdgesPerVertex,Nfile);

  Nmax = 2.^SCALE;                           # Max vertex ID.
  M = EdgesPerVertex .* Nmax;                # Total number of edges.
  myFiles = collect(1:Nfile).';              # Set list of files.
  #
  tab = Char(9)
  nl = Char(10)
  Niter = 20                                      # Number of PageRank iterations.
  c = 0.15                                        # PageRank damping factor.

  println("Number of Edges: " * string(M) * ", Maximum Possible Vertex: " * string(Nmax));


  ########################################################
  # Kernel 0: Generate a Graph500 Kronecker graph and save to data files.
  ########################################################
  println("Kernel 0: Generate Graph, Write Edges");
  tic();
  let SCALE=SCALE,EdgesPerVertex=EdgesPerVertex,Nfile=Nfile
    pmap(i -> begin
      fname = "data/K0/" * string(i) * ".tsv";
      println("  Writing: " * fname);                          # Read filename.
      srand(i);                                                # Set random seed to be unique for this file.
      ut, vt = KronGraph500NoPerm(SCALE,EdgesPerVertex./Nfile);  # Generate data.
      writeuv(fname, ut, vt)
      nothing
    end, myFiles)
  end
  K0time = toq();
  println("K0 Time: " * string(K0time) * ", Edges/sec: " * string(M./K0time));

  ########################################################
  # Kernel 1: Read data, sort data, and save to files.
  ########################################################
  println("Kernel 1: Read, Sort, Write Edges");
  tic();

  # Each worker sorts a range of Nmax. Local range computed from a small sample.
    u=SharedArray(Int, M)
    v=SharedArray(Int, M)

    # Read in all the files into one array.
    nPerFile = div(M,Nfile)
    pmap(i->begin
      fname = "data/K0/" * string(i) * ".tsv"
      println("  Reading: " * fname);  # Read filename.
      ut,vt = StrFileRead(fname)
      # Concatenate to u,v
      startOffset = (i-1)*nPerFile + 1
      endOffSet = i*nPerFile

      u[startOffset:endOffSet] = ut
      v[startOffset:endOffSet] = vt
      nothing
    end, myFiles)

  K1time1 = toq();
  tic();

    # Take a small sample and calculate how to distribute the range among workers.
    sample_size = M>4096 ? 4096 : M
    sample = u[1:sample_size]
    sort!(sample)
    lower_bounds = sample[[1+(x-1)*div(sample_size, Nfile) for x in 1:Nfile]]
    upper_bounds = [x-1 for x in lower_bounds[2:end]]
    push!(upper_bounds, M)

   localsort(i, lb, ub) = begin
      idxs = find(x->(x >= lb && x <= ub), u)   # To be sorted on this worker
      sortset_local = u[idxs]
      sorted_local_idxs = sortperm(sortset_local)
      sorted_u_idxs = Int[idxs[x] for x in sorted_local_idxs]

      # write out both u and v
      fname = "data/K1/" * string(i) * ".tsv"
      println("  Writing: " * fname)                              # Create filename.
      writeuv(fname, u[sorted_u_idxs], v[sorted_u_idxs])
      (string(i) * ".tsv", length(sorted_u_idxs))
    end

    results=pmap(localsort, myFiles, lower_bounds, upper_bounds)

    # write out metadata
    f=open("data/K1/METADATA", "w")
    for (fn, l) in results
      println(f, fn, ",", l)
    end
    close(f)


  K1time2 = toq();
  K1time = K1time1 + K1time2;
  println("K1 Time (reading):" * string(K1time1) * ", Edges/sec: " * string(M./K1time1));
  println("K1 Time (sorting and writing):" * string(K1time2) * ", Edges/sec: " * string(M./K1time2));
  println("K1 Time: " * string(K1time) * ", Edges/sec: " * string(M./K1time));

# Serial sort impl.
#    sortIndex = sortperm(u)                      # Sort starting vertices.
#    u = u[sortIndex]                                  # Get starting vertices.
#    v = v[sortIndex]                                  # Get ending vertices.

  return K0time,K1time,0.0,0.0

  ########################################################
  # Kernel 2: Read data, filter data.
  ########################################################
  println("Kernel 2: Read, Filter Edges");
  tic();
    # Read in all the files into one array.
    for i in myFiles
      fname = "data/K1/" * string(i) * ".tsv";
      println("  Reading: " * fname);                # Read filename.
      ut,vt = StrFileRead(fname);
      if i == 1
         u = ut; v = vt;                             # Initialize starting and ending vertices.
      else
         append!(u, ut)
         append!(v, vt)
         # Get the rest of starting and ending vertices.
      end
    end

    # Construct adjacency matrix.
    A = sparse(vec(u),vec(v),1.0,Nmax,Nmax)      # Create adjacency matrix.

    # Filter and weight the adjacency matrix.
    din = sum(A,1)                               # Compute in degree.
    A[find(din == maximum(din))]=0               # Eliminate the super-node.
    A[find(din == 1)]=0                          # Eliminate the leaf-node.
    dout = sum(A,2)                              # Compute the out degree.
    is = find(dout)                               # Find vertices with outgoing edges (dout > 0).
    DoutInvD = zeros(Nmax)        # Create diagonal weight matrix.
    DoutInvD[is] = 1./dout[is]
    scale!(DoutInvD, A)           # Apply weight matrix.
  K2time = toq();
  println("K2 Time: " * string(K2time) * ", Edges/sec: " * string(M./K2time));


  ########################################################
  # Kernel 3: Compute PageRank.
  ########################################################
  println("Kernel 3: PageRank");
  tic();

    r = rand(1,Nmax);                     # Generate a random starting rank.
    r = r ./ norm(r,1);                   # Normalize
    a = (1-c) ./ Nmax;                    # Create damping vector

    for i=1:Niter
        s = r * A
        scale!(s, c)
        r = s .+ (a * sum(r,2));                # Compute PageRank.
    end

    r=r./norm(r,1);

  K3time = toq();
  println("  Sum of PageRank: " * string(sum(r)) );     # Force all computations to occur.
  println("K3 Time: " * string(K3time) * ", Edges/sec: " * string(Niter.*M./K3time));

  return K0time,K1time,K2time,K3time

end

########################################################
# PageRank Pipeline Benchmark
# Architect: Dr. Jeremy Kepner (kepner@ll.mit.edu)
# Julia Translation: Dr. Chansup Byun (cbyun@ll.mit.edu)
# MIT
########################################################
# (c) <2015> Massachusetts Institute of Technology
########################################################
