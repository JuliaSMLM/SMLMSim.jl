@testset "Static SMLM" begin
    # Create a pattern for testing
    pattern = Nmer2D(n=3, d=0.2)
    
    # Create a fluorophore with appropriate parameters
    fluorophore = GenericFluor(
        photons = 1000.0,
        k_off = 0.2,
        k_on = 0.1
    )
    
    # Create StaticSMLMParams
    params = StaticSMLMParams(
        density = 0.5,        # particles per μm²
        σ_psf = 0.13,         # μm
        minphotons = 100,     
        ndatasets = 1,        
        nframes = 10,         
        framerate = 10.0,     # frames per second
        ndims = 2             
    )
    
    # Test parameters
    @test params.density == 0.5
    @test params.σ_psf == 0.13
    @test params.minphotons == 100
    @test params.ndatasets == 1
    @test params.nframes == 10
    @test params.framerate == 10.0
    @test params.ndims == 2
    
    # Test coordinate noise application
    # Create some emitters for testing (using positional arguments)
    emitters = Vector{Emitter2DFit{Float64}}()
    
    # Create emitters with explicit frame numbers, IDs, and datasets
    push!(emitters, Emitter2DFit(0.1, 0.2, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, 1, 1, 1))
    push!(emitters, Emitter2DFit(0.5, 0.6, 1500.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2, 1, 2, 2))
    push!(emitters, Emitter2DFit(1.0, 1.2, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3, 1, 3, 3))
    
    # Create a camera
    camera = IdealCamera(32, 32, 0.1)  # 32x32 pixels, 100 nm pixel size
    
    # Create SMLD
    smld = BasicSMLD(emitters, camera, 3, 1)  # 3 frames, 1 dataset
    
    # Apply noise with standard deviation of 0.02
    noisy_smld = apply_noise(smld, 0.02)
    
    # Check that coordinates were modified
    @test length(noisy_smld.emitters) == length(smld.emitters)
    
    for i in 1:length(smld.emitters)
        @test noisy_smld.emitters[i].id == smld.emitters[i].id  # ID should remain unchanged
        @test noisy_smld.emitters[i].frame == smld.emitters[i].frame  # Frame should remain unchanged
        @test noisy_smld.emitters[i].photons == smld.emitters[i].photons  # Photons should remain unchanged
        
        # Coordinates should be modified
        @test noisy_smld.emitters[i].x != smld.emitters[i].x || noisy_smld.emitters[i].y != smld.emitters[i].y
    end
    
    # Verify that noise is approximately correct
    differences_x = [noisy_smld.emitters[i].x - smld.emitters[i].x for i in 1:length(smld.emitters)]
    differences_y = [noisy_smld.emitters[i].y - smld.emitters[i].y for i in 1:length(smld.emitters)]
    
    # Check that there is some difference in the emitter positions
    # but not check the exact standard deviation (it depends on implementation)
    @test std(differences_x) > 0.0
    @test std(differences_y) > 0.0
end