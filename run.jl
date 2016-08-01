#!/usr/bin/env julia

using ArgParse
using BenchmarkTools

const kernel_args = Dict{Int, Vector{Any}}(
    0 => Any["/tmp/kernel0.tsv", 20, 16]
)

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

        mod = @eval PageRank.$(Symbol(target))
        state = @eval $mod.setup()

        # Process all kernels
        for k in 0:3
            isdefined(mod, Symbol("kernel$k")) || break
            kernel = @eval $mod.$(Symbol("kernel$(k)"))
            setup = @eval $mod.$(Symbol("kernel$(k)_setup"))
            teardown = @eval $mod.$(Symbol("kernel$(k)_teardown"))

            if length(parsed_args["kernel"]) == 0 || in(parsed_args["kernel"], k)
                println("Benchmarking kernel $k...")
                benchmark = @benchmarkable $kernel($state, args) setup=(
                    args = $setup($state, $(kernel_args[0]...))
                ) teardown=(
                    $teardown($state, $args)
                )
                trials = run(benchmark)

                if parsed_args["verbose"]
                    println(trials)
                else
                    println("Executed in ", BenchmarkTools.prettytime(time(trials)))
                end
            else
                println("Running kernel $k...")
                @eval begin
                    args = $setup($state, $(kernel_args[0]...))
                    $kernel($state, args)
                    @eval $teardown($state, $args)
                end
            end
        end

        @eval $mod.teardown($state)
        println()
    end
end

main(ARGS)
