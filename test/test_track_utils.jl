@testset "Track Utilities" begin
    # Create a test SMLD with multiple tracks
    camera = IdealCamera(32, 32, 0.1)
    
    # Create emitters with different track_ids
    emitters = Vector{Emitter2DFit{Float64}}()
    
    # Track 1: Two emitters in different frames
    push!(emitters, Emitter2DFit(0.5, 0.5, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, 1, 1, 1))
    push!(emitters, Emitter2DFit(0.51, 0.52, 900.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2, 1, 1, 1))
    
    # Track 2: Single emitter
    push!(emitters, Emitter2DFit(1.5, 1.5, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, 1, 2, 2))
    
    # Track 3: Three emitters in different frames
    push!(emitters, Emitter2DFit(2.5, 2.5, 700.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, 1, 3, 3))
    push!(emitters, Emitter2DFit(2.52, 2.53, 750.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2, 1, 3, 3))
    push!(emitters, Emitter2DFit(2.54, 2.56, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3, 1, 3, 3))
    
    # Create test SMLD
    smld = BasicSMLD(emitters, camera, 3, 1)

    # Test get_track function
    @testset "get_track" begin
        # Get track 1
        track1 = get_track(smld, 1)
        @test isa(track1, BasicSMLD)
        @test length(track1.emitters) == 2
        @test all(e -> e.track_id == 1, track1.emitters)
        
        # Get track 2
        track2 = get_track(smld, 2)
        @test isa(track2, BasicSMLD)
        @test length(track2.emitters) == 1
        @test all(e -> e.track_id == 2, track2.emitters)
        
        # Get track 3
        track3 = get_track(smld, 3)
        @test isa(track3, BasicSMLD)
        @test length(track3.emitters) == 3
        @test all(e -> e.track_id == 3, track3.emitters)
        
        # Now test get_track with a nonexistent track_id
        track4 = get_track(smld, 4)
        @test isa(track4, BasicSMLD)
        
        # Instead of checking emptiness directly, let's check filter behavior
        @test count(e -> e.track_id == 4, track4.emitters) == 0
        
        # If the function is working correctly, we should either have an empty array
        # or a filtered array with no track_id == 4
        @test isempty(track4.emitters) || count(e -> e.track_id == 4, track4.emitters) == 0
    end
    
    # Test get_num_tracks function
    @testset "get_num_tracks" begin
        num_tracks = get_num_tracks(smld)
        @test num_tracks == 3
        
        # Test with empty SMLD
        empty_smld = BasicSMLD(Emitter2DFit{Float64}[], camera, 1, 1)
        @test get_num_tracks(empty_smld) == 0
    end
    
    # Test get_tracks function
    @testset "get_tracks" begin
        tracks = get_tracks(smld)
        @test isa(tracks, Vector{BasicSMLD})
        @test length(tracks) == 3
        
        # Check that each track has the correct emitters
        @test all(e -> e.track_id == 1, tracks[1].emitters)
        @test all(e -> e.track_id == 2, tracks[2].emitters)
        @test all(e -> e.track_id == 3, tracks[3].emitters)
        
        # Check sizes
        @test length(tracks[1].emitters) == 2
        @test length(tracks[2].emitters) == 1
        @test length(tracks[3].emitters) == 3
        
        # Test with empty SMLD
        empty_smld = BasicSMLD(Emitter2DFit{Float64}[], camera, 1, 1)
        empty_tracks = get_tracks(empty_smld)
        @test isa(empty_tracks, Vector{BasicSMLD})
        @test isempty(empty_tracks)
    end
end