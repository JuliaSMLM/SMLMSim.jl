using SMLMSim
using Test
using Distributions
using LinearAlgebra
using SMLMData
using Statistics
using MicroscopePSFs

@testset "Diffusion SMLM" begin
    # Create diffusion simulation parameters
    params = DiffusionSMLMParams(
        density = 0.5,             # molecules per μm²
        box_size = 5.0,            # 5μm box
        diff_monomer = 0.1,        # μm²/s
        diff_dimer = 0.05,         # μm²/s 
        diff_dimer_rot = 0.1,      # rad²/s
        k_off = 1.0,               # s⁻¹
        r_react = 0.015,           # reaction radius in μm
        d_dimer = 0.03,            # monomer separation in dimer in μm
        dt = 0.001,                # time step in s
        t_max = 0.1,               # short simulation for testing
        boundary = "periodic",     # periodic boundary conditions
        ndims = 2,                 # 2D simulation
        camera_framerate = 100.0,  # frames per second
        camera_exposure = 0.01     # exposure time in s
    )
    
    # Test parameters
    @test params.density == 0.5
    @test params.box_size == 5.0
    @test params.diff_monomer == 0.1
    @test params.diff_dimer == 0.05
    @test params.diff_dimer_rot == 0.1
    @test params.k_off == 1.0
    @test params.r_react == 0.015
    @test params.d_dimer == 0.03
    @test params.dt == 0.001
    @test params.t_max == 0.1
    @test params.boundary == "periodic"
    @test params.ndims == 2
    @test params.camera_framerate == 100.0
    @test params.camera_exposure == 0.01
    
    # Test simulation (using a smaller system for faster tests)
    small_params = DiffusionSMLMParams(
        density = 0.5,             # molecules per μm²
        box_size = 2.0,            # 2μm box for faster tests
        diff_monomer = 0.1,        # μm²/s
        diff_dimer = 0.05,         # μm²/s
        diff_dimer_rot = 0.1,      # rad²/s
        k_off = 1.0,               # s⁻¹
        r_react = 0.015,           # μm
        d_dimer = 0.03,            # μm
        dt = 0.001,                # s
        t_max = 0.1,               # s
        boundary = "periodic",     # periodic boundary conditions
        ndims = 2,                 # 2D simulation
        camera_framerate = 100.0,  # frames per second
        camera_exposure = 0.01     # s
    )
    
    # Calculate expected number of molecules based on density and box size
    expected_n_molecules = round(Int, small_params.density * small_params.box_size^2)
    
    # Run simulation
    result = simulate(small_params)
    
    # From our earlier check, we know result is a BasicSMLD, not a structure with trajectories
    # Let's test the SMLD properties instead
    
    # Test that we have emitters in our simulation result
    @test !isempty(result)
    
    # Test metadata
    @test haskey(result.metadata, "simulation_type")
    
    # Check that emitters have proper physical units
    if !isempty(result.emitters)
        e = result.emitters[1]
        @test e.x >= 0 && e.x <= small_params.box_size
        @test e.y >= 0 && e.y <= small_params.box_size
        if small_params.ndims == 3
            @test e.z >= 0 && e.z <= small_params.box_size
        end
    end
    
    # Test camera integration with diffusion results
    if !isempty(result) && result.camera !== nothing
        # Create a PSF model
        psf = GaussianPSF(0.13)  # 130 nm PSF width
        
        # Generate camera images
        images = gen_images(result, psf)
        
        # Verify we have at least one frame
        @test size(images, 3) > 0
    end
end