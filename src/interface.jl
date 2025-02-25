"""
    simulate(;
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

# Examples
```julia
# Basic 2D simulation with default parameters
camera = IdealCamera(1:128, 1:128, 0.1)  # 100nm pixels
smld_true, smld_model, smld_noisy = simulate(camera=camera)

# 3D simulation with custom parameters
camera = IdealCamera(1:256, 1:256, 0.1)
pattern = Nmer3D(n=6, d=0.2)
smld_true, smld_model, smld_noisy = simulate(
    ρ=2.0,                # 2 patterns per μm²
    σ_psf=0.15,           # 150nm PSF width
    minphotons=100,       # minimum photons for detection
    ndatasets=5,          # 5 independent datasets
    nframes=2000,         # 2000 frames per dataset
    framerate=100.0,      # 100 frames per second
    pattern=pattern,
    molecule=GenericFluor(; γ=1e5, q=[0 20; 5e-2 0]),
    camera=camera,
    zrange=[-2.0, 2.0]    # 4μm axial range
)
```

# Note
For 3D patterns, PSF width is automatically adjusted for axial dimension according
to typical microscope parameters.
"""
function simulate(;
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
    # Input validation
    if ρ <= 0
        throw(ArgumentError("Density (ρ) must be positive"))
    end
    
    if σ_psf <= 0
        throw(ArgumentError("PSF width (σ_psf) must be positive"))
    end
    
    if minphotons < 0
        throw(ArgumentError("Minimum photon count must be non-negative"))
    end
    
    if ndatasets <= 0
        throw(ArgumentError("Number of datasets must be positive"))
    end
    
    if nframes <= 0
        throw(ArgumentError("Number of frames must be positive"))
    end
    
    if framerate <= 0
        throw(ArgumentError("Frame rate must be positive"))
    end
    
    if length(zrange) != 2 || zrange[1] >= zrange[2]
        throw(ArgumentError("zrange must be a vector of two values [min_z, max_z] where min_z < max_z"))
    end
    
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
        "pattern_type" => string(typeof(pattern)),
        "pattern_params" => Dict(
            fn => getfield(pattern, fn) 
            for fn in fieldnames(typeof(pattern)) 
            if fn ∉ [:x, :y, :z]
        )
    )
    smld_true = BasicSMLD(emitters, camera, 1, 1, metadata)

    # Apply kinetic model
    smld_model = kinetic_model(smld_true, molecule, nframes, framerate;
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

"""
    simulate(pattern::Pattern; kwargs...)

Simulate SMLM data using a specific pattern.

This is a convenience method that passes the pattern directly to the main simulate function.
All other parameters are specified as keyword arguments.

# Arguments
- `pattern::Pattern`: The pattern to use for simulation

# Keyword Arguments
Same as the main `simulate` function

# Returns
- `Tuple{BasicSMLD, BasicSMLD, BasicSMLD}`: (true_positions, model_kinetics, noisy_data)

# Example
```julia
pattern = Line2D(λ=5.0, endpoints=[(-2.0, 0.0), (2.0, 0.0)])
smld_true, smld_model, smld_noisy = simulate(pattern, ρ=0.5, nframes=2000)
```
"""
function simulate(pattern::Pattern; kwargs...)
    simulate(; pattern=pattern, kwargs...)
end

# Maintain backward compatibility with the old sim function
"""
    sim(; kwargs...)

Alias for `simulate` function. Provided for backward compatibility.

See `simulate` for full documentation and parameters.
"""
function sim(; kwargs...)
    @warn "The `sim` function is deprecated, use `simulate` instead."
    simulate(; kwargs...)
end

function sim(pattern::Pattern; kwargs...)
    @warn "The `sim` function is deprecated, use `simulate` instead."
    simulate(pattern; kwargs...)
end