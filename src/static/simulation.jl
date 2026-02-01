"""
    Main implementation for static SMLM simulations
"""

"""
    simulate(params::StaticSMLMParams;
             starting_conditions::Union{Nothing, SMLD, Vector{<:AbstractEmitter}}=nothing,
             pattern::Pattern=nothing,
             labeling::AbstractLabeling=FixedLabeling(),
             molecule::Molecule=GenericFluor(photons=1e4, k_off=50.0, k_on=1e-2),
             camera::AbstractCamera=IdealCamera(1:128, 1:128, 0.1),
             state1::Union{Int, Symbol}=:equilibrium,
             burn_in::Real=0.0)

Generate simulated static SMLM data with realistic blinking kinetics
and localization uncertainty.

# Arguments
- `params::StaticSMLMParams`: Simulation parameters
- `starting_conditions::Union{Nothing, SMLD, Vector{<:AbstractEmitter}}`: Optional starting conditions instead of generating patterns
- `pattern::Pattern`: Pattern to use (default depends on params.ndims)
- `labeling::AbstractLabeling`: Labeling strategy for fluorophore attachment (default: FixedLabeling() = 1 per site)
- `molecule::Molecule`: Fluorophore model for blinking simulation
- `camera::AbstractCamera`: Camera model for detection simulation
- `state1::Union{Int, Symbol}=:equilibrium`: Initial fluorophore state:
  - `::Int`: Specific state (1=on, 2=off typically)
  - `:equilibrium`: Sample from equilibrium distribution (default)
  For models with photobleaching (absorbing states), use `state1=1`.
- `burn_in::Real=0.0`: Pre-illumination time in seconds before recording starts.
  Simulates experimental protocol where laser is on before data collection,
  allowing some molecules to bleach and reach pseudo-equilibrium.

# Returns
- `Tuple{BasicSMLD, BasicSMLD, BasicSMLD}`: (true_positions, model_kinetics, noisy_data)
    - true_positions: Ground truth emitter positions (after labeling)
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

# Run with Poisson labeling (average 1.5 fluorophores per binding site)
smld_true, smld_model, smld_noisy = simulate(params;
    pattern=pattern,
    labeling=PoissonLabeling(1.5)
)

# Run with partial labeling efficiency (80% of sites labeled)
smld_true, smld_model, smld_noisy = simulate(params;
    pattern=pattern,
    labeling=PoissonLabeling(1.0; efficiency=0.8)
)

# Run with photobleaching model and 5s burn-in
k_off, k_on, k_bleach = 10.0, 1.0, 0.1
Q = [-(k_off+k_bleach) k_off k_bleach; k_on -k_on 0.0; 0.0 0.0 0.0]
fluor = GenericFluor(1e4, Q)
smld_true, smld_model, smld_noisy = simulate(params;
    molecule=fluor,
    state1=1,        # Start in ON state (required for absorbing states)
    burn_in=5.0      # 5 seconds pre-illumination
)

# Run with custom starting conditions (labeling not applied)
custom_emitters = [
    Emitter2DFit{Float64}(x, y, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; σ_xy=0.0, track_id=i)
    for (i, (x, y)) in enumerate(zip(rand(10), rand(10)))
]
smld_true, smld_model, smld_noisy = simulate(params; starting_conditions=custom_emitters)
```
# Note
- The `params.σ_psf` value is used directly for lateral uncertainty (σx, σy) in both 2D and 3D.
- For 3D simulations, the axial uncertainty (σz) is scaled by a factor of 3 (i.e., σz = 3 * σ_psf).
- If `starting_conditions` is provided, it will be used instead of generating patterns, and labeling is not applied.
- Labeling and molecule are separate concepts: labeling controls how many fluorophores per binding site,
  molecule controls the photophysics (blinking) of each fluorophore.
"""
function simulate(params::StaticSMLMParams;
                 starting_conditions::Union{Nothing, SMLD, Vector{<:AbstractEmitter}}=nothing,
                 pattern::Union{Pattern,Nothing}=nothing,
                 labeling::AbstractLabeling=FixedLabeling(),
                 molecule::Molecule=GenericFluor(photons=1e4, k_off=50.0, k_on=1e-2),
                 camera::AbstractCamera=IdealCamera(1:128, 1:128, 0.1),
                 state1::Union{Int, Symbol}=:equilibrium,
                 burn_in::Real=0.0)
    
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

        # Generate binding site coordinates based on pattern type
        # uniform2D/uniform3D now return pattern_ids as well
        coords_with_ids = if pattern isa Pattern2D
            x, y, pattern_ids = uniform2D(params.density, pattern, field_x, field_y)
            ((x, y), pattern_ids)
        else
            x, y, z, pattern_ids = uniform3D(params.density, pattern, field_x, field_y; zrange=params.zrange)
            ((x, y, z), pattern_ids)
        end
        coords, pattern_ids = coords_with_ids

        # Apply labeling to expand binding sites to fluorophore positions
        coords, pattern_ids = apply_labeling(coords, pattern_ids, labeling)

        # Create emitters for true positions
        emitters = if pattern isa Pattern2D
            x, y = coords
            [Emitter2DFit{Float64}(
                x[i], y[i],          # positions
                1.0,                 # nominal photons
                0.0,                 # background
                0.0, 0.0,            # σ_x, σ_y
                0.0, 0.0;            # σ_photons, σ_bg
                σ_xy=0.0,            # x-y covariance (0 for symmetric PSF)
                frame=1,
                dataset=1,
                track_id=i,
                id=pattern_ids[i]
            ) for i in 1:length(x)]
        else
            x, y, z = coords
            [Emitter3DFit{Float64}(
                x[i], y[i], z[i],    # positions
                1.0,                 # nominal photons
                0.0,                 # background
                0.0, 0.0, 0.0,       # σ_x, σ_y, σ_z
                0.0, 0.0;            # σ_photons, σ_bg
                σ_xy=0.0,            # x-y covariance (0 for symmetric PSF)
                frame=1,
                dataset=1,
                track_id=i,
                id=pattern_ids[i]
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
        metadata["labeling_type"] = string(typeof(labeling))
        metadata["labeling_params"] = Dict(
            fn => getfield(labeling, fn)
            for fn in fieldnames(typeof(labeling))
        )
        
        # Create SMLD with true positions
        smld_true = BasicSMLD(emitters, camera, 1, 1, metadata)
    end

    # Apply kinetic model
    smld_model = kinetic_model(smld_true, molecule, params.nframes, params.framerate;
                             ndatasets=params.ndatasets, minphotons=params.minphotons,
                             state1=state1, burn_in=burn_in)

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
