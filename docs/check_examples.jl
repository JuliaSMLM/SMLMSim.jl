using SMLMSim
using Documenter
using CairoMakie
using MicroscopePSFs
using Distributions

function check_example_file(filename)
    println("Testing examples in $filename...")
    md_content = read(filename, String)
    
    # Extract all doctests
    doctest_blocks = []
    for m in eachmatch(r"```jldoctest.*?```"s, md_content)
        push!(doctest_blocks, m.match)
    end
    
    # Print information
    println("Found $(length(doctest_blocks)) doctest blocks")
    
    for (i, block) in enumerate(doctest_blocks)
        println("\nTesting block $i:")
        
        # Extract the code (remove backticks and jldoctest header)
        code = replace(block, r"```jldoctest.*?\n"s => "")
        code = replace(code, r"```\s*$" => "")
        code = replace(code, r"#\s*output.*?\n" => "")
        
        try
            # Evaluate the code in the current module
            include_string(Main, code)
            println("✓ Block $i passed")
        catch e
            println("✗ Block $i failed with error:")
            println(e)
            return false
        end
    end
    
    return true
end

# Test static examples
println("\n== Testing Static Examples ==")
check_example_file(joinpath(@__DIR__, "src/static/examples.md"))

# Test diffusion examples
println("\n== Testing Diffusion Examples ==")
check_example_file(joinpath(@__DIR__, "src/diffusion/examples.md"))