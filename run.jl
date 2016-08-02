#!/usr/bin/env julia

using ArgParse
using BenchmarkTools


#=

Implementation interface
------------------------

- make sure $NAME/pagerank.jl exists and defines module PageRank$(ucfirst(NAME))
- provide setup() returning a global state, and teardown() to destroy it
  if you don't need global state, just return `nothing`
- provide kernel$n(...) methods with following function signature requirements:
  - last argument represents the global state
    if you don't care about state, provide a default value of `nothing` for ease of REPL use
  - kernel0(output_directory, scale, edge_factor, scale=see_previous_point)
  - other kernels: return values of previous kernel is passed as-is (along with the state)

=#

function main(args)
    s = ArgParseSettings(description = "PageRank benchmark driver")
    @add_arg_table s begin
        "--kernel","-k"
            action = :append_arg
            arg_type = Int
            help = "kernels to benchmark (defaults to all kernels)"
        "--verbose","-v"
            action = :store_true
            help = "verbose output"
        "implementations"
            nargs = '+'
            required = true
            help = "implementations to benchmark"
    end
    parsed_args = parse_args(s)

    for impl in parsed_args["implementations"]
        target = String(impl)
        println("Benchmarking $target\n")

        entrypoint = joinpath(Base.source_dir(), target, "pagerank.jl")
        isfile(entrypoint) || error("Implementation $impl does not have an entrypoint pagerank.jl")
        include(entrypoint)

        mod = @eval $(Symbol("PageRank$(ucfirst(target))"))
        state = @eval $mod.setup()

        # Process all kernels
        kernel_args = ("/tmp", 20, 16)
        for k in 0:3
            isdefined(mod, Symbol("kernel$k")) || break
            kernel = @eval $mod.$(Symbol("kernel$(k)"))

            # Run the kernel
            println("Running kernel $k...")
            kernel_expr = :( $kernel($(kernel_args...), $state) )
            kernel_output = @eval $kernel_expr

            # Benchmark this kernel if requested
            if length(parsed_args["kernel"]) == 0 || in(parsed_args["kernel"], k)
                print("Benchmarking kernel $k... ")
                benchmark = @benchmarkable $kernel_expr
                trials = run(benchmark)

                if parsed_args["verbose"]
                    println(trials)
                else
                    println(BenchmarkTools.prettytime(time(trials)))
                end
            end

            kernel_args = kernel_output
        end

        @eval $mod.teardown($state)
        println()
    end
end

main(ARGS)
