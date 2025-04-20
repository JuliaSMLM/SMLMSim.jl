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
    # Test 2D pattern types
    @testset "2D Patterns" begin
        # Test Nmer2D
        pattern2d = Nmer2D(n=3, d=0.2)  # d is the diameter of the pattern
        
        # Verify we have 3 positions
        @test pattern2d.n == 3
        @test length(pattern2d.x) == 3
        @test length(pattern2d.y) == 3
        
        # Test that all points are at approximately the same distance from the origin
        distances = sqrt.(pattern2d.x.^2 .+ pattern2d.y.^2)
        @test all(isapprox.(distances, 0.1, atol=1e-10))  # d/2 = 0.1
        
        # Test Line2D
        line2d = Line2D(λ=10.0, endpoints=[(-0.25, 0.0), (0.25, 0.0)])
        @test length(line2d.x) == line2d.n
        @test length(line2d.y) == line2d.n
        @test line2d.λ == 10.0  # Linear density
        
        # Test that points are within the line endpoints
        @test all(x -> -0.25 <= x <= 0.25, line2d.x)
        @test all(isapprox.(line2d.y, 0.0, atol=1e-10))
        
        # Test pattern rotation
        rotated_pattern = Nmer2D(n=3, d=0.2)
        original_x = copy(rotated_pattern.x)
        original_y = copy(rotated_pattern.y)
        
        # Rotate 90 degrees - explicitly qualify the rotate! function
        SMLMSim.rotate!(rotated_pattern, π/2)
        
        # Check that coordinates have changed
        @test !all(isapprox.(rotated_pattern.x, original_x, atol=1e-10))
        @test !all(isapprox.(rotated_pattern.y, original_y, atol=1e-10))
        
        # Check that the distance from origin is preserved (rotation should preserve distances)
        original_distances = sqrt.(original_x.^2 .+ original_y.^2)
        rotated_distances = sqrt.(rotated_pattern.x.^2 .+ rotated_pattern.y.^2)
        @test all(isapprox.(original_distances, rotated_distances, atol=1e-10))
    end
    
    # Test 3D pattern types
    @testset "3D Patterns" begin
        # Test Nmer3D
        pattern3d = Nmer3D(n=4, d=0.3)  # d is the diameter of the pattern
        
        # Verify we have 4 positions
        @test pattern3d.n == 4
        @test length(pattern3d.x) == 4
        @test length(pattern3d.y) == 4
        @test length(pattern3d.z) == 4
        
        # Test that points form a circle in xy plane (z=0)
        distances_xy = sqrt.(pattern3d.x.^2 .+ pattern3d.y.^2)
        @test all(isapprox.(distances_xy, 0.3/2, atol=1e-10))  # d/2 = 0.15
        @test all(isapprox.(pattern3d.z, 0.0, atol=1e-10))  # All z values should be 0
        
        # Test Line3D
        line3d = Line3D(λ=8.0, endpoints=[(-0.3, 0.0, -0.2), (0.3, 0.0, 0.2)])
        @test length(line3d.x) == line3d.n
        @test length(line3d.y) == line3d.n
        @test length(line3d.z) == line3d.n
        @test line3d.λ == 8.0  # Linear density
        
        # Test that points are between the endpoints
        @test all(x -> -0.3 <= x <= 0.3, line3d.x)
        @test all(isapprox.(line3d.y, 0.0, atol=1e-10))
        @test all(z -> -0.2 <= z <= 0.2, line3d.z)
        
        # Test 3D rotation with Euler angles
        rotated_3d = Nmer3D(n=4, d=0.3)
        original_x = copy(rotated_3d.x)
        original_y = copy(rotated_3d.y)
        original_z = copy(rotated_3d.z)
        
        # Rotate using Euler angles 
        SMLMSim.rotate!(rotated_3d, π/4, π/6, π/3)  # 45°, 30°, 60° rotations
        
        # Check that coordinates have changed
        @test !all(isapprox.(rotated_3d.x, original_x, atol=1e-10))
        @test !all(isapprox.(rotated_3d.y, original_y, atol=1e-10))
        @test !all(isapprox.(rotated_3d.z, 0.0, atol=1e-10))  # Z should no longer be 0
        
        # Check that the distance from origin is preserved
        original_distances = sqrt.(original_x.^2 .+ original_y.^2 .+ original_z.^2)
        rotated_distances = sqrt.(rotated_3d.x.^2 .+ rotated_3d.y.^2 .+ rotated_3d.z.^2)
        @test all(isapprox.(original_distances, rotated_distances, atol=1e-10))
    end
    
    # Test pattern distribution functions
    @testset "Pattern Distribution" begin
        # Test uniform2D
        nmer = Nmer2D(n=3, d=0.2)
        x, y = SMLMSim.uniform2D(0.5, nmer, 10.0, 10.0)
        @test length(x) == length(y)
        @test length(x) % 3 == 0  # Should be multiple of 3 (pattern size)
        # Note: Coordinates may extend slightly beyond field boundaries due to pattern placement
        # rather than strict bounds testing, check they're reasonably close
        @test all(xi -> -1.0 <= xi <= 11.0, x)  # Allow modest extension beyond field
        @test all(yi -> -1.0 <= yi <= 11.0, y)  # Allow modest extension beyond field
        
        # Test uniform3D 
        nmer3d = Nmer3D(n=4, d=0.3)
        x, y, z = SMLMSim.uniform3D(0.5, nmer3d, 10.0, 10.0, zrange=[-1.0, 1.0])
        @test length(x) == length(y) == length(z)
        @test length(x) % 4 == 0  # Should be multiple of 4 (pattern size)
        # Similar relaxed bounds for 3D
        @test all(xi -> -1.0 <= xi <= 11.0, x)
        @test all(yi -> -1.0 <= yi <= 11.0, y)
        @test all(zi -> -1.5 <= zi <= 1.5, z)  # Allow modest extension beyond z range
    end
end