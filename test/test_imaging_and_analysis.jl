using SMLMSim.InteractionDiffusion
using MicroscopePSFs
using Random

function test_imaging_and_analysis()
    @testset "Image Generation" begin
        # Create a simple system with one molecule
        camera = IdealCamera(1:50, 1:50, 0.1)  # 5μm x 5μm field with 100nm pixels
        emitter = Emitter2D{Float64}(2.5, 2.5, 1000.0)  # Center of field
        molecule = DiffusingMolecule(emitter, 1, 1, nothing, false)
        
        system = DiffusingMoleculeSystem(
            [molecule], camera, 5.0, 1, 1, Dict{String,Any}()
        )
        
        # Create PSF and generate image
        psf = GaussianPSF(0.15)
        img = gen_image(psf, system, 1, photons=1000.0, bg=5.0, poisson_noise=false)
        
        # Test basic image properties
        @test size(img) == (length(camera.pixel_edges_y)-1, length(camera.pixel_edges_x)-1)
        @test all(img .>= 5.0)  # At least background everywhere
        
        # Test image sequence generation
        # Create a minimal simulation
        params = DiffusionSMLMParams(
            density = 1.0,
            box_size = 2.0,
            dt = 0.01,
            t_max = 0.03  # Just 3 frames
        )
        systems = simulate(params)
        
        # Generate image sequence
        imgs = gen_image_sequence(psf, systems, 
                                 photons=1000.0, bg=5.0, 
                                 frame_integration=1, poisson_noise=false)
        
        # Test basic sequence properties
        @test size(imgs, 3) == length(systems)
        
        # Test frame integration
        imgs_integrated = gen_image_sequence(psf, systems, 
                                            photons=1000.0, bg=5.0, 
                                            frame_integration=3, poisson_noise=false)
        
        @test size(imgs_integrated, 3) == 1  # Should have just one frame with integration=3
    end
    
    @testset "Dimer Imaging" begin
        # Create a system with dimers
        params = DiffusionSMLMParams(
            density = 5.0,
            box_size = 1.0,
            dt = 0.01,
            t_max = 0.01  # Just one frame
        )
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
        
        # Create a sequence with just this one system
        systems = [system]
        
        # Create PSF and generate dimer images
        psf = GaussianPSF(0.15)
        dimer_images = gen_dimer_images(systems, psf, 
                                       photons=1000.0, bg=5.0,
                                       frame_integration=1, poisson_noise=false)
        
        # Test basic image properties
        @test size(dimer_images, 3) == 1
    end
    
    @testset "Full Simulation Pipeline" begin
        # Test the complete workflow with minimal simulation
        # This is a crucial test as it ensures all components work together
        
        # Create simulation parameters
        params = DiffusionSMLMParams(
            density = 1.0,
            box_size = 1.0,
            dt = 0.01,
            t_max = 0.02  # Just 2 frames
        )
        
        # Create PSF
        psf = GaussianPSF(0.15)
        
        # Run simulation and imaging
        images, systems = simulate_and_image(params, psf,
            photons=1000.0, bg=5.0,
            frame_integration=1, poisson_noise=false)
        
        # Test basic outputs
        @test isa(images, Array{Float64, 3})
        @test length(systems) == 2
        @test size(images, 3) == 2
    end
    
    @testset "Core Interface" begin
        # Test the main simulate interfaces with StaticSMLMParams
        camera = IdealCamera(1:20, 1:20, 0.1)
        
        # Basic interface with StaticSMLMParams
        params = StaticSMLMParams(
            density=1.0,
            nframes=5  # Minimal
        )
        smld_true, smld_model, smld_noisy = simulate(params, camera=camera)
        
        # Test basic outputs
        @test smld_true isa BasicSMLD
        @test smld_model isa BasicSMLD
        @test smld_noisy isa BasicSMLD
        
        # With explicit pattern
        pattern = Nmer2D(n=4, d=0.1)
        params = StaticSMLMParams(
            density=1.0,
            nframes=5
        )
        smld_true, smld_model, smld_noisy = simulate(params, pattern=pattern, camera=camera)
        
        @test smld_true isa BasicSMLD
    end
    
    # Skip visualization tests in normal test runs as they're slow and mostly check for errors
    if get(ENV, "RUN_VISUALIZATION_TESTS", "false") == "true"
        @testset "Visualization" begin
            # Create a small simulation
            params = DiffusionSMLMParams(
                density = 2.0,
                box_size = 1.0,
                dt = 0.01,
                t_max = 0.03  # Just 3 frames
            )
            systems = simulate(params)
            
            # Test with a temporary file (no assertion, just checking for errors)
            temp_file = "test_vis_$(rand(1:1000)).mp4"
            
            try
                visualize_sequence(systems, filename=temp_file, framerate=10)
            finally
                if isfile(temp_file)
                    rm(temp_file)
                end
            end
        end
    end
end