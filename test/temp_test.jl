using SMLMSim
using Test

@testset "Diffusion helpers refactoring" begin
    # Test that we can still use types
    camera = IdealCamera(1:100, 1:100, 0.1)
    emitter = Emitter2D(1.0, 2.0, 1000.0)
    mol1 = DiffusingMolecule(emitter, 1, 1, nothing, false)
    
    # Create a second molecule
    emitter2 = Emitter2D(1.1, 2.1, 1000.0)
    mol2 = DiffusingMolecule(emitter2, 1, 2, nothing, false)
    
    # Create a system with both molecules
    system = DiffusingMoleculeSystem(
        [mol1, mol2], 
        camera, 
        10.0,  # box size
        1,     # frames
        1,     # datasets
        Dict{String,Any}("simulation_type" => "diffusion")
    )
    
    # Test helper functions
    r = calc_r(mol1, mol2)
    @test isa(r, Float64)
    @test r ≈ sqrt(0.1^2 + 0.1^2)
    
    ϕ = calc_ϕ(mol1, mol2)
    @test isa(ϕ, Float64)
    
    θ = calc_θ(mol1, mol2)
    @test isa(θ, Float64)
    @test θ ≈ π/2
    
    # Test dimerize!
    dimerize!(mol1, mol2, 0.05)
    @test mol1.state == 2
    @test mol2.state == 2
    @test mol1.link == 2
    @test mol2.link == 1
    
    # Calculate new distance after dimerization
    r_after = calc_r(mol1, mol2)
    @test r_after ≈ 0.05
    
    # Test monomerize!
    monomerize!(mol1, system)
    @test mol1.state == 1
    @test mol2.state == 1
    @test isnothing(mol1.link)
    @test isnothing(mol2.link)
end

println("All tests passed!")
