using SMLMSim
using Test
using Distributions
using LinearAlgebra

# Main test set
@testset "SMLMSim.jl" begin
    # Include individual test files
    @testset "Core Module" begin
        include("test_core.jl")
    end
    
    @testset "Static SMLM" begin
        include("test_static.jl")
    end
    
    @testset "Diffusion SMLM" begin
        include("test_diffusion.jl")
    end
    
    @testset "Camera Images" begin
        include("test_camera_images.jl")
    end
    
    @testset "Integration Tests" begin
        include("test_integration.jl")
    end
end