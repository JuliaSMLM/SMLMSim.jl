using SMLMSim
using Test

@testset "Analysis module refactoring" begin
    # Create a test system
    camera = IdealCamera(1:100, 1:100, 0.1)
    
    # Create molecules, some as monomers, some as dimers
    molecules = []
    
    # Create 8 molecules total, 4 monomers and 2 dimers
    for i in 1:8
        emitter = Emitter2D(rand(), rand(), 1000.0)
        molecule = DiffusingMolecule(emitter, 1, i, nothing, false)
        push!(molecules, molecule)
    end
    
    # Make two dimers (4 molecules)
    molecules[1].state = 2
    molecules[1].link = 2
    molecules[2].state = 2
    molecules[2].link = 1
    
    molecules[3].state = 2
    molecules[3].link = 4
    molecules[4].state = 2
    molecules[4].link = 3
    
    # Create system
    system = DiffusingMoleculeSystem(
        molecules, 
        camera, 
        10.0,  # box size
        1,     # frames
        1,     # datasets
        Dict{String,Any}("simulation_type" => "diffusion")
    )
    
    # Create a sequence of systems
    systems = [system for _ in 1:5]
    
    # Test get_dimers for a single system
    dimer_system = get_dimers(system)
    @test length(dimer_system.molecules) == 4
    @test all(mol -> mol.state == 2, dimer_system.molecules)
    
    # Test get_dimers for a sequence
    dimer_systems = get_dimers(systems)
    @test length(dimer_systems) == 5
    @test all(sys -> length(sys.molecules) == 4, dimer_systems)
    
    # Test analyze_dimer_fraction
    fractions = analyze_dimer_fraction(systems)
    @test length(fractions) == 5
    @test all(f -> f ≈ 0.5, fractions)  # 4/8 = 0.5
    
    # Test analyze_dimer_lifetime
    # This just checks that the function doesn't error, as it's not fully implemented
    @test_logs (:warn, "Dimer lifetime analysis not yet fully implemented") analyze_dimer_lifetime(systems)
end

println("All analysis tests passed!")
