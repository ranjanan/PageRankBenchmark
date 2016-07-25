using BufferedStreams

"""
Optimized function that writes list of edges as tab seperated values
Tries to be perform as little allocation as possible
"""
function writetsv(file, ij1, ij2)
   arr = Vector{UInt8}(64)
   for (i, j) in zip(ij1, ij2)
      nchars = dec!(arr, i)
      writeto(file, arr, nchars)
      write(file, UInt8('\t'))
      nchars = dec!(arr, j)
      writeto(file, arr, nchars)
      write(file, UInt8('\n'))
   end
end

function writeto(io, arr, nchars)
   @inbounds for i in 1:nchars
      write(io, arr[i])
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

# Assumes positive number
function dec!(arr::Vector{UInt8}, x::Int64)
   nchars = i = Base.ndigits0z(x)
   if i > length(arr) # only resize if array is to small
      resize!(arr, i)
   end
   while i > 0
      arr[i] = '0'+rem(x,10)
      x ÷= 10
      i -= 1
   end
   nchars
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

