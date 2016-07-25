using BufferedStreams

"""
Optimized function that writes list of edges as tab seperated values
Tries to be perform as little allocation as possible
"""
function writetsv(file, ij1, ij2)
   for (i, j) in zip(ij1, ij2)
      write(file, dec(i))
      write(file, '\t')
      write(file, dec(j))
      write(file, '\n')
   end
end

function readtsv(file)
   ij1 = Array{Int64}(0)
   ij2 = Array{Int64}(0)
   const TAB = UInt8('\t')
   const NL = UInt8('\n')
   while !eof(file)
      push!(ij1, parseuntil(file,TAB))
      push!(ij2, parseuntil(file,NL))
   end
   ij1, ij2
end

# Exploit the domain knowledge we have
# Only positive numbers
# Taken from CSV.jl
function parseuntil(file, delim::UInt8)
   v = 0
   b = read(file, UInt8)
   const ZERO = UInt8('0')
   while b != delim
      v *= 10
      v += b - ZERO
      b = read(file, UInt8)
   end
   v
end

