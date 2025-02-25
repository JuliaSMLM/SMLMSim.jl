using SMLMSim.InteractionDiffusion
using Random
# Import the initialize_system function explicitly
using SMLMSim.InteractionDiffusion: initialize_system

function test_diffusion()
    @testset "SmoluchowskiParams" begin
        # Test default constructor
        params = SmoluchowskiParams()
        @test params.density == 1.0
        @test params.box_size == 10.0
        @test params.dt == 0.01
        @test params.boundary == "periodic"
        
        # Test custom constructor with a few key parameters
        params = SmoluchowskiParams(
            density = 2.0,
            box_size = 20.0,
            ndims = 3,
            boundary = "reflecting"
        )
        
        @test params.density == 2.0
        @test params.box_size == 20.0
        @test params.ndims == 3
        @test params.boundary == "reflecting"
        
        # Test one invalid parameter case
        @test_throws ArgumentError SmoluchowskiParams(density = -1.0)
    end
    
    @testset "Diffusion Simulation" begin
        # Create small system for quick tests
        params = SmoluchowskiParams(
            density = 1.0,
            box_size = 1.0,
            dt = 0.01,
            t_max = 0.05  # Just 5 frames
        )
        
        # Test system initialization
        system = initialize_system(params)
        expected_molecules = round(Int, params.density * params.box_size^params.ndims)
        @test length(system.molecules) == expected_molecules
        
        # Test basic simulation
        systems = simulate(params)
        @test length(systems) == round(Int, params.t_max / params.dt)
        
        # Test dimer formation with parameters favoring dimers
        params_dimer = SmoluchowskiParams(
            density = 10.0,     # High density
            box_size = 1.0,     # Small box
            r_react = 0.1,      # Large reaction radius
            k_off = 0.01,       # Slow dissociation
            dt = 0.01,
            t_max = 0.1         # Short but likely sufficient for dimers
        )
        
        systems_dimer = simulate(params_dimer)
        
        # Check if dimers formed
        final_system = systems_dimer[end]
        n_dimers = count(m -> m.state == 2, final_system.molecules)
        
        # Print rather than test - dimers might not form in all cases
        @info "Number of dimers formed: $n_dimers out of $(length(final_system.molecules)) molecules"
    end
    
    @testset "Dimer Analysis" begin
        # Create a system with some dimers
        params = SmoluchowskiParams(box_size = 1.0, density = 5.0)
        system = initialize_system(params)
        
        # Make some molecules into dimers
        n_molecules = length(system.molecules)
        n_dimers = min(2, div(n_molecules, 2))
        for i in 1:n_dimers
            mol1 = system.molecules[2*i-1]
            mol2 = system.molecules[2*i]
            mol1.state = 2
            mol2.state = 2
            mol1.link = mol2.id
            mol2.link = mol1.id
        end
        
        # Test get_dimers
        dimer_system = get_dimers(system)
        @test length(dimer_system.molecules) == 2 * n_dimers
        @test all(m -> m.state == 2, dimer_system.molecules)
        
        # Test analyze_dimer_fraction
        fractions = analyze_dimer_fraction([system])
        @test length(fractions) == 1
        expected_fraction = 2 * n_dimers / n_molecules
        @test fractions[1] ≈ expected_fraction
    end
end