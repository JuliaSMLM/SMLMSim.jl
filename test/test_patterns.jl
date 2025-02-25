function test_patterns()
    @testset "Pattern Types" begin
        # Test Nmer2D as representative 2D pattern
        nmer = Nmer2D(n=6, d=0.2)
        @test nmer.n == 6
        @test nmer.d == 0.2
        @test length(nmer.x) == 6
        
        # Test points are on a circle of specified diameter (just check a few points)
        @test sqrt(nmer.x[1]^2 + nmer.y[1]^2) ≈ nmer.d/2 atol=1e-10
        
        # Basic existence test for Line2D
        line2d = Line2D(λ=5.0)
        @test line2d.λ == 5.0
        @test length(line2d.x) > 0
        
        # Test Nmer3D as representative 3D pattern
        nmer3d = Nmer3D(n=8, d=0.5)
        @test nmer3d.n == 8 
        @test length(nmer3d.z) == 8
        @test all(nmer3d.z .== 0.0)  # Should all be in z=0 plane
        
        # Basic existence test for Line3D
        line3d = Line3D()
        @test length(line3d.z) > 0
    end
    
    @testset "Pattern Generation and Rotation" begin
        # Test uniform2D with basic checks
        nmer = Nmer2D(n=4, d=0.1)
        
        # Calculate the maximum possible extent of the pattern
        # For Nmer2D, this is half the diameter
        pattern_radius = nmer.d/2
        
        # Field dimensions
        field_x = 10.0
        field_y = 10.0
        
        x, y = uniform2D(2.0, nmer, field_x, field_y)
        
        # Basic checks
        @test length(x) == length(y)
        @test length(x) % nmer.n == 0
        
        # Check ranges with pattern extent considered
        # Pattern centers are in [0, field_x], points can be at most pattern_radius away
        @test all(-pattern_radius .<= x .<= field_x + pattern_radius)
        @test all(-pattern_radius .<= y .<= field_y + pattern_radius)
        
        # Test rotation changes coordinates
        x_orig = copy(nmer.x)
        y_orig = copy(nmer.y)
        
        # Use fully qualified name for rotate!
        SMLMSim.rotate!(nmer, π/2)
        
        # Just check that something changed - no need for exact calculations
        @test any(nmer.x .!= x_orig) || any(nmer.y .!= y_orig)
        
        # Basic test for 3D functions
        nmer3d = Nmer3D(n=4, d=0.2)
        pattern_radius_3d = nmer3d.d/2
        field_x_3d = 5.0
        field_y_3d = 5.0
        zrange = [-1.0, 1.0]  # Default z-range
        
        x, y, z = uniform3D(1.0, nmer3d, field_x_3d, field_y_3d)
        @test length(x) == length(y) == length(z)
        
        # Check ranges with pattern extent considered
        @test all(-pattern_radius_3d .<= x .<= field_x_3d + pattern_radius_3d)
        @test all(-pattern_radius_3d .<= y .<= field_y_3d + pattern_radius_3d)
        @test all(zrange[1] - pattern_radius_3d .<= z .<= zrange[2] + pattern_radius_3d)
    end
end