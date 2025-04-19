@testset "Core - CTMC" begin
    # Test CTMC state transitions
    # Rate matrix: state 1 -> state 2 at rate 2.0, state 2 -> state 1 at rate 1.0
    q = [-2.0 2.0; 1.0 -1.0]  # Proper rate matrix (diagonals are negative, rows sum to zero)
    simulation_time = 10.0
    initial_state = 1
    
    ctmc = CTMC(q, simulation_time, initial_state)
    
    # Test initial state
    @test get_state(ctmc, 0.0) == initial_state
    
    # Create a CTMC with known transitions for deterministic testing
    custom_times = [0.0, 1.0, 3.0, 6.0]
    custom_states = [1, 2, 1, 2]
    custom_ctmc = CTMC(simulation_time, custom_times, custom_states)
    
    # Test states at various time points
    @test get_state(custom_ctmc, 0.5) == 1  # Between 0 and 1, should be state 1
    @test get_state(custom_ctmc, 1.0) == 2  # At transition point, should be state 2
    @test get_state(custom_ctmc, 2.0) == 2  # Between 1 and 3, should be state 2
    @test get_state(custom_ctmc, 5.0) == 1  # Between 3 and 6, should be state 1
    @test get_state(custom_ctmc, 9.0) == 2  # After last transition, should be state 2
    
    # Test get_next function
    next_state, next_time = get_next(custom_ctmc, 0.5)
    @test next_state == 2
    @test next_time == 1.0
    
    next_state, next_time = get_next(custom_ctmc, 4.0)
    @test next_state == 2
    @test next_time == 6.0
end

@testset "Core - Molecules" begin
    # Test GenericFluor creation and properties with current constructor
    # Keywords are: photons, k_off, k_on
    fluor = GenericFluor(
        photons = 1000.0,
        k_off = 0.2,
        k_on = 0.1
    )
    
    @test fluor.γ == 1000.0
    @test size(fluor.q) == (2, 2)
    @test fluor.q[1, 2] == 0.2  # k_off
    @test fluor.q[2, 1] == 0.1  # k_on
    
    # Test fluorophore properties directly
    @test isa(fluor.q, Matrix)
    @test size(fluor.q) == (2, 2)  # 2x2 matrix for a simple fluorophore
end

@testset "Core - Patterns" begin
    # Test 2D pattern creation
    pattern2d = Nmer2D(n=3, d=0.2)  # d is the diameter of the pattern
    
    # Verify we have 3 positions
    @test pattern2d.n == 3
    @test length(pattern2d.x) == 3
    @test length(pattern2d.y) == 3
    
    # Test that all points are at approximately the same distance from the origin
    distances = sqrt.(pattern2d.x.^2 .+ pattern2d.y.^2)
    @test all(isapprox.(distances, 0.1, atol=1e-10))  # d/2 = 0.1
    
    # Test Line2D if it exists
    if isdefined(SMLMSim, :Line2D)
        line = Line2D(λ=10.0, endpoints=[(-0.25, 0.0), (0.25, 0.0)])
        @test length(line.x) == line.n
        @test length(line.y) == line.n
    end
    
    # Test uniform2D if it exists
    if isdefined(SMLMSim, :uniform2D)
        nmer = Nmer2D(n=3, d=0.2)
        x, y = uniform2D(0.5, nmer, 10.0, 10.0)
        @test length(x) == length(y)
        # We can't reliably test bounds as this depends on implementation details
    end
end