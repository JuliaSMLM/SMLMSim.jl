"""
    Main implementation for static SMLM simulations
"""

"""
    simulate(params::StaticSMLMParams; 
             starting_conditions::Union{Nothing, SMLD, Vector{<:AbstractEmitter}}=nothing,
             pattern::Pattern=nothing,
             molecule::Molecule=GenericFluor(photons=1e4, k_off=50.0, k_on=1e-2),
             camera::AbstractCamera=IdealCamera(1:128, 1:128, 0.1))

Generate simulated static SMLM data with realistic blinking kinetics
and localization uncertainty.

# Arguments
- `params::StaticSMLMParams`: Simulation parameters
- `starting_conditions::Union{Nothing, SMLD, Vector{<:AbstractEmitter}}`: Optional starting conditions instead of generating patterns
- `pattern::Pattern`: Pattern to use (default depends on params.ndims)
- `molecule::Molecule`: Fluorophore model for blinking simulation
- `camera::AbstractCamera`: Camera model for detection simulation

# Returns
- `Tuple{BasicSMLD, BasicSMLD, BasicSMLD}`: (true_positions, model_kinetics, noisy_data)
    - true_positions: Ground truth emitter positions
    - model_kinetics: Positions with simulated blinking
    - noisy_data: Positions with blinking and localization uncertainty

# Example
```julia
# Create parameters
params = StaticSMLMParams(
    density = 2.0,              # 2 patterns per μm²
    σ_psf = 0.15,         # 150nm PSF width
    minphotons = 100,     # 100 photons for detection
    ndatasets = 5,        # 5 independent datasets
    nframes = 2000,       # 2000 frames
    framerate = 100.0     # 100 frames per second
)

# Run simulation with Nmer pattern
pattern = Nmer3D(n=6, d=0.2)
smld_true, smld_model, smld_noisy = simulate(params; pattern=pattern)

# Run with custom starting conditions
custom_emitters = [
    Emitter2DFit{Float64}(x, y, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; track_id=i)
    for (i, (x, y)) in enumerate(zip(rand(10), rand(10)))
]
smld_true, smld_model, smld_noisy = simulate(params; starting_conditions=custom_emitters)
```
# Note
- The `params.σ_psf` value is used directly for lateral uncertainty (σx, σy) in both 2D and 3D.
- For 3D simulations, the axial uncertainty (σz) is scaled by a factor of 3 (i.e., σz = 3 * σ_psf).
- If `starting_conditions` is provided, it will be used instead of generating patterns.
"""
function simulate(params::StaticSMLMParams; 
                 starting_conditions::Union{Nothing, SMLD, Vector{<:AbstractEmitter}}=nothing,
                 pattern::Union{Pattern,Nothing}=nothing,
                 molecule::Molecule=GenericFluor(photons=1e4, k_off=50.0, k_on=1e-2),
                 camera::AbstractCamera=IdealCamera(1:128, 1:128, 0.1))
    
    # Initialize metadata
    metadata = Dict{String,Any}(
        "simulation_parameters" => params
    )
    
    # Process starting conditions if provided
    if starting_conditions !== nothing
        # Extract emitters from starting_conditions
        emitters = if starting_conditions isa SMLD
            # Deep copy the emitters from the SMLD
            deepcopy.(starting_conditions.emitters)
        else
            # Already a vector of emitters
            deepcopy.(starting_conditions)
        end
        
        # Validate that emitters have appropriate type
        emitter_type = eltype(emitters)
        if !(emitter_type <: AbstractEmitter)
            error("Starting conditions must contain valid emitters")
        end
        
        # Set up metadata for ground truth with starting conditions
        metadata["simulation_type"] = "ground_truth_from_starting_conditions"
        metadata["source"] = "user_provided"
        
        # Create SMLD with true positions from starting conditions
        smld_true = BasicSMLD(emitters, camera, 1, 1, metadata)
    else
        # Use pattern-based generation (original code path)
        
        # Use appropriate default pattern if none provided
        if pattern === nothing
            pattern = params.ndims == 3 ? Nmer3D() : Nmer2D()
        end
        
        # Get field size in microns from camera
        centers_x, centers_y = get_pixel_centers(camera)
        field_x = maximum(centers_x) - minimum(centers_x)
        field_y = maximum(centers_y) - minimum(centers_y)

        # Generate coordinates based on pattern type
        coords = if pattern isa Pattern2D
            uniform2D(params.density, pattern, field_x, field_y)
        else
            uniform3D(params.density, pattern, field_x, field_y; zrange=params.zrange)
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

        # Add pattern-specific metadata
        metadata["simulation_type"] = "ground_truth"
        metadata["density"] = params.density
        metadata["pattern_type"] = string(typeof(pattern))
        metadata["pattern_params"] = Dict(
            fn => getfield(pattern, fn) 
            for fn in fieldnames(typeof(pattern)) 
            if fn ∉ [:x, :y, :z]
        )
        
        # Create SMLD with true positions
        smld_true = BasicSMLD(emitters, camera, 1, 1, metadata)
    end

    # Apply kinetic model
    smld_model = kinetic_model(smld_true, molecule, params.nframes, params.framerate;
                             ndatasets=params.ndatasets, minphotons=params.minphotons)

    # Determine PSF scaling based on emitter dimensionality
    emitter_type = eltype(smld_true.emitters)
    σ_scaled = if emitter_type <: Emitter2DFit
        params.σ_psf  # already in microns
    else
        [params.σ_psf, params.σ_psf, 3.0*params.σ_psf]  # typical 3D PSF scaling
    end
    
    # Add localization noise with appropriate PSF scaling
    smld_noisy = apply_noise(smld_model, σ_scaled)

    return smld_true, smld_model, smld_noisy
end
