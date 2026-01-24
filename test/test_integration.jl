@testset "Integration Tests" begin
    @testset "Static SMLM Workflow" begin
        # Create StaticSMLMParams
        params = StaticSMLMParams(
            density = 0.5,        # particles per μm²
            σ_psf = 0.13,         # μm
            minphotons = 100,     
            ndatasets = 1,        
            nframes = 5,          # small number for testing
            framerate = 10.0,     # frames per second
            ndims = 2             
        )
        
        # Create some emitters for testing (using Emitter2DFit)
        emitters = Vector{Emitter2DFit{Float64}}()
        
        # Create emitters with explicit frame numbers
        push!(emitters, Emitter2DFit{Float64}(0.5, 0.5, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0; σ_xy=0.0, frame=1, dataset=1, track_id=1, id=1))
        push!(emitters, Emitter2DFit{Float64}(1.5, 1.5, 2000.0, 0.0, 0.0, 0.0, 0.0, 0.0; σ_xy=0.0, frame=1, dataset=1, track_id=2, id=2))
        push!(emitters, Emitter2DFit{Float64}(0.5, 0.5, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0; σ_xy=0.0, frame=2, dataset=1, track_id=3, id=3))
        
        # Create a camera for testing
        camera = IdealCamera(32, 32, 0.1)  # 32x32 pixels, 100 nm pixel size
        
        # Create SMLD container
        smld = BasicSMLD(emitters, camera, 2, 1)  # 2 frames, 1 dataset
        
        # Add localization uncertainty with standard deviation of 0.02
        noisy_smld = apply_noise(smld, 0.02)
        
        # Create a PSF model
        psf = GaussianPSF(0.13)  # 130 nm PSF width
        
        # Generate camera images
        images = gen_images(noisy_smld, psf; bg=10.0)
        
        # Verify workflow outputs
        @test length(noisy_smld.emitters) > 0
        @test size(images, 3) == 2  # Should have 2 frames
    end
    
    @testset "Diffusion SMLM Workflow" begin
        # Create diffusion simulation parameters (small system for quick tests)
        params = DiffusionSMLMParams(
            density = 0.5,             # molecules per μm²
            box_size = 2.0,            # 2μm box
            diff_monomer = 0.1,        # μm²/s
            diff_dimer = 0.05,         # μm²/s
            diff_dimer_rot = 0.1,      # rad²/s
            k_off = 1.0,               # s⁻¹
            r_react = 0.015,           # μm
            d_dimer = 0.03,            # μm
            dt = 0.001,                # s
            t_max = 0.1,               # short simulation for testing
            boundary = "periodic",     # periodic boundary conditions
            ndims = 2,                 # 2D simulation
            camera_framerate = 100.0,  # frames per second
            camera_exposure = 0.01     # s
        )
        
        # Run simulation
        smld = simulate(params)
        
        # Test that we have emitters
        @test !isempty(smld.emitters)
        
        # Create a PSF model if we have emitters and a camera
        if !isempty(smld.emitters) && smld.camera !== nothing
            # Create a PSF model
            psf = GaussianPSF(0.13)  # 130 nm PSF width
            
            # Generate camera images
            images = gen_images(smld, psf; bg=10.0)
            
            # Verify we have images
            @test size(images, 3) > 0
        end
        
        # Verify outputs
        @test smld.n_frames > 0
        @test haskey(smld.metadata, "simulation_type")
    end
end