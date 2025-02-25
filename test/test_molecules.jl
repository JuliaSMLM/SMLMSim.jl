using SMLMSim.InteractionDiffusion: DiffusingMolecule, DiffusingMoleculeSystem

function test_molecules()
    @testset "GenericFluor" begin
        # Test various constructors with minimal checks
        # Use Float64 values for q to match AbstractFloat requirement
        fluor_default = GenericFluor(γ=1e5, q=[1.0])  # Use explicit Float64 value
        @test fluor_default.γ == 1e5
        
        # Custom parameters
        fluor_custom = GenericFluor(2e4, [0.0 5.0; 0.1 0.0])  # Use explicit Float64 values
        @test fluor_custom.γ == 2e4
        @test fluor_custom.q[1,2] == 5.0
        
        # Keyword constructor
        fluor_kw = GenericFluor(γ=5e4, q=[0.0 10.0; 0.2 0.0])  # Use explicit Float64 values
        @test fluor_kw.γ == 5e4
    end
    
    @testset "DiffusingMolecule" begin
        # Create a 2D emitter and molecule
        emitter = Emitter2D{Float64}(1.0, 2.0, 1000.0)
        molecule = DiffusingMolecule(emitter, 1, 1, nothing, false)
        
        # Test basic properties
        @test molecule.state == 1
        @test molecule.id == 1
        @test molecule.link === nothing
        
        # Test property access and modification
        @test molecule.x == 1.0
        molecule.x = 3.0
        @test molecule.x == 3.0
        @test molecule.emitter.x == 3.0  # Should change underlying emitter too
    end
    
    @testset "DiffusingMoleculeSystem" begin
        # Create minimal system
        camera = IdealCamera(1:100, 1:100, 0.1)
        emitter = Emitter2D{Float64}(1.0, 2.0, 1000.0)
        molecule = DiffusingMolecule(emitter, 1, 1, nothing, false)
        
        # Create system
        system = DiffusingMoleculeSystem(
            [molecule], camera, 10.0, 1, 1, Dict{String,Any}("test" => "value")
        )
        
        # Test basic properties
        @test length(system.molecules) == 1
        @test system.box_size == 10.0
        @test length(system.emitters) == 1
        @test system.emitters[1] === molecule.emitter
    end
end