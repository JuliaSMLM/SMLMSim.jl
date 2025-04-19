using SMLMSim
using Test
using SMLMData
using Distributions
using LinearAlgebra
using MicroscopePSFs

# Include test files
include("test_patterns.jl")
include("test_molecules.jl")
include("test_kinetics.jl")
include("test_diffusion.jl")
include("test_imaging_and_analysis.jl")

@testset "SMLMSim.jl" begin
    @testset "Patterns" begin
        test_patterns()
    end

    @testset "Molecules" begin
        test_molecules()
    end

    @testset "Kinetic Modeling" begin
        test_kinetics()
    end

    @testset "Diffusion" begin
        test_diffusion()
    end

    @testset "Imaging & Analysis" begin
        test_imaging_and_analysis()
    end
end