"""
    sim(;
        ρ=1.0,                   # density in particles per square micron
        σ_psf=0.13,              # PSF width in microns
        minphotons=50,           # minimum photons for detection
        ndatasets=10,            # number of datasets to simulate
        nframes=1000,            # number of frames per dataset
        framerate=50.0,          # frames per second
        pattern::Pattern=Nmer2D(),
        molecule::Molecule=GenericFluor(; q=[0 50; 1e-2 0]),
        camera::AbstractCamera=IdealCamera(1:128, 1:128, 0.1),  # 100nm pixels
        zrange::Vector{Float64}=[-1.0, 1.0]  # axial range in microns
    )

Generate simulated single molecule localization data with realistic blinking kinetics
and localization uncertainty.

# Returns
- `Tuple{BasicSMLD, BasicSMLD, BasicSMLD}`: (true_positions, model_kinetics, noisy_data)
    - true_positions: Ground truth emitter positions
    - model_kinetics: Positions with simulated blinking
    - noisy_data: Positions with blinking and localization uncertainty

# Details
1. Generates emitter positions based on pattern and density
2. Simulates blinking using kinetic model
3. Adds localization uncertainty based on photon counts
4. All coordinates are in physical units (microns)

# Note
For 3D patterns, PSF width is automatically adjusted for axial dimension according
to typical microscope parameters.
"""
function sim(;
    ρ=1.0,
    σ_psf=0.13,
    minphotons=50,
    ndatasets=10,
    nframes=1000,
    framerate=50.0,
    pattern::Pattern=Nmer2D(),
    molecule::Molecule=GenericFluor(; q=[0 50; 1e-2 0]),
    camera::AbstractCamera=IdealCamera(1:128, 1:128, 0.1),
    zrange::Vector{Float64}=[-1.0, 1.0] 
)
    # Get field size in microns from camera
    centers_x, centers_y = get_pixel_centers(camera)
    field_x = maximum(centers_x) - minimum(centers_x)
    field_y = maximum(centers_y) - minimum(centers_y)

    # Generate coordinates based on pattern type
    coords = if pattern isa Pattern2D
        uniform2D(ρ, pattern, field_x, field_y)
    else
        uniform3D(ρ, pattern, field_x, field_y; zrange=zrange)
    end

    # Create emitters for true positions
    emitters = if pattern isa Pattern2D
        x, y = coords
        [Emitter2DFit{Float64}(
            x[i], y[i],          # positions
            1.0,                 # nominal photons
            0.0,                 # background
            0.0, 0.0,            # σ_x, σ_y
            0.0, 0.0;            # σ_photons, σ_bg
            frame=1,
            dataset=1,
            track_id=i
        ) for i in 1:length(x)]
    else
        x, y, z = coords
        [Emitter3DFit{Float64}(
            x[i], y[i], z[i],    # positions
            1.0,                 # nominal photons
            0.0,                 # background
            0.0, 0.0, 0.0,       # σ_x, σ_y, σ_z
            0.0, 0.0;            # σ_photons, σ_bg
            frame=1,
            dataset=1,
            track_id=i
        ) for i in 1:length(x)]
    end

    # Create SMLD with true positions
    metadata = Dict{String,Any}(
        "simulation_type" => "ground_truth",
        "density" => ρ,
        "pattern" => string(typeof(pattern))
    )
    smld_true = BasicSMLD(emitters, camera, 1, 1, metadata)

    # Apply kinetic model
    smld_model = kineticmodel(smld_true, molecule, nframes, framerate;
                             ndatasets=ndatasets, minphotons=minphotons)

    # Add localization noise with appropriate PSF scaling
    σ_scaled = if pattern isa Pattern2D
        σ_psf  # already in microns
    else
        [σ_psf, σ_psf, 3.0*σ_psf]  # typical 3D PSF scaling
    end
    
    smld_noisy = noise(smld_model, σ_scaled)

    return smld_true, smld_model, smld_noisy
end