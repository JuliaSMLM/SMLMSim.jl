"""
    Main implementation for static SMLM simulations
"""

"""
    simulate(params::StaticSMLMConfig;
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
- `params::StaticSMLMConfig`: Simulation parameters
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
- `Tuple{BasicSMLD, SimInfo}`: (noisy_data, info)
    - noisy_data: Positions with blinking and localization uncertainty
    - info: SimInfo containing smld_true, smld_model, timing, and counts

# Example
```julia
# Create parameters
params = StaticSMLMConfig(
    density = 2.0,              # 2 patterns per μm²
    σ_psf = 0.15,         # 150nm PSF width
    minphotons = 100,     # 100 photons for detection
    ndatasets = 5,        # 5 independent datasets
    nframes = 2000,       # 2000 frames
    framerate = 100.0     # 100 frames per second
)

# Run simulation with Nmer pattern
pattern = Nmer3D(n=6, d=0.2)
smld_noisy, info = simulate(params; pattern=pattern)

# Access intermediate results
smld_true = info.smld_true
smld_model = info.smld_model

# Run with Poisson labeling (average 1.5 fluorophores per binding site)
smld_noisy, info = simulate(params;
    pattern=pattern,
    labeling=PoissonLabeling(1.5)
)

# Run with partial labeling efficiency (80% of sites labeled)
smld_noisy, info = simulate(params;
    pattern=pattern,
    labeling=PoissonLabeling(1.0; efficiency=0.8)
)

# Run with photobleaching model and 5s burn-in
k_off, k_on, k_bleach = 10.0, 1.0, 0.1
Q = [-(k_off+k_bleach) k_off k_bleach; k_on -k_on 0.0; 0.0 0.0 0.0]
fluor = GenericFluor(1e4, Q)
smld_noisy, info = simulate(params;
    molecule=fluor,
    state1=1,        # Start in ON state (required for absorbing states)
    burn_in=5.0      # 5 seconds pre-illumination
)

# Run with custom starting conditions (labeling not applied)
custom_emitters = [
    Emitter2DFit{Float64}(x, y, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; σ_xy=0.0, track_id=i)
    for (i, (x, y)) in enumerate(zip(rand(10), rand(10)))
]
smld_noisy, info = simulate(params; starting_conditions=custom_emitters)
```
# Note
- The field is the camera's full extent, from its first to its last pixel edge in x and in y.
  The sample continues past it: pattern centers are placed uniformly, at `params.density`
  patterns per μm², over the field widened on every side by the pattern's radius r (the largest
  distance of a binding site from the pattern's center), and only binding sites that land on
  the camera are kept. The emitter density on the camera is therefore flat to its edges, and a
  pattern of zero extent (`Nmer2D(n=1, d=0.0)`, `Nmer3D(n=1, d=0.0)`) is placed over the field
  itself. `info.n_patterns` counts the patterns with a binding site on the camera, before
  labeling (a pattern whose sites all get zero fluorophores still counts), and their ids run
  1, 2, ....
- Through 0.7.2 the field was the span of the pixel centers, starting at 0: one pixel short in
  x and in y, so the last pixel row and column held no pattern center; and binding sites of a
  pattern centered near an edge were kept off the camera, with none coming in from outside it,
  so the density was lower within r of each edge. Seeded runs differ from 0.7.2.
- The `params.σ_psf` value is used directly for lateral uncertainty (σx, σy) in both 2D and 3D.
- For 3D simulations, the axial uncertainty (σz) is scaled by a factor of 3 (i.e., σz = 3 * σ_psf).
- If `starting_conditions` is provided, it will be used instead of generating patterns, and labeling is not applied.
- Labeling and molecule are separate concepts: labeling controls how many fluorophores per binding site,
  molecule controls the photophysics (blinking) of each fluorophore.
"""
function simulate(params::StaticSMLMConfig;
                 starting_conditions::Union{Nothing, SMLD, Vector{<:AbstractEmitter}}=nothing,
                 pattern::Union{Pattern,Nothing}=nothing,
                 labeling::AbstractLabeling=FixedLabeling(),
                 molecule::Molecule=GenericFluor(photons=1e4, k_off=50.0, k_on=1e-2),
                 camera::AbstractCamera=IdealCamera(1:128, 1:128, 0.1),
                 state1::Union{Int, Symbol}=:equilibrium,
                 burn_in::Real=0.0)

    start_time = time_ns()
    n_patterns = 0

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
        n_patterns = 0  # User-provided, no patterns generated
    else
        # Use pattern-based generation (original code path)

        # Use appropriate default pattern if none provided
        if pattern === nothing
            pattern = params.ndims == 3 ? Nmer3D() : Nmer2D()
        end

        # The camera's full extent, from its first to its last pixel edge (SMLMData's
        # pixel-edge convention)
        x_lo, x_hi = Float64(first(camera.pixel_edges_x)), Float64(last(camera.pixel_edges_x))
        y_lo, y_hi = Float64(first(camera.pixel_edges_y)), Float64(last(camera.pixel_edges_y))

        # The sample continues past the camera: pattern centers are drawn over the extent
        # widened on every side by the pattern's radius r, so every center that can put a
        # binding site on the camera is drawn, and only sites that land on the camera are kept.
        # The emitter density is then flat to the edge. A pattern of zero extent (r = 0) is
        # drawn over the extent itself and keeps every site.
        r = _pattern_radius(pattern)
        x_origin, y_origin = x_lo - r, y_lo - r
        field_x = (x_hi - x_lo) + 2r
        field_y = (y_hi - y_lo) + 2r

        # Generate binding site coordinates based on pattern type, with pattern centers in
        # [0, field], then shift them to the widened field's origin
        coords_with_ids = if pattern isa Pattern2D
            x, y, pattern_ids = uniform2D(params.density, pattern, field_x, field_y)
            ((x, y), pattern_ids)
        else
            x, y, z, pattern_ids = uniform3D(params.density, pattern, field_x, field_y; zrange=params.zrange)
            ((x, y, z), pattern_ids)
        end
        coords, pattern_ids = coords_with_ids
        coords[1] .+= x_origin
        coords[2] .+= y_origin

        if r > 0
            # Keep the binding sites on the camera, and number the patterns that kept one 1, 2, ...
            keep = [x_lo <= coords[1][i] <= x_hi && y_lo <= coords[2][i] <= y_hi
                    for i in eachindex(pattern_ids)]
            coords = map(c -> c[keep], coords)
            renumber = Dict{Int,Int}()
            pattern_ids = [get!(renumber, id, length(renumber) + 1) for id in pattern_ids[keep]]
        end

        # Count patterns (those with a binding site on the camera) before labeling expands them
        n_patterns = length(unique(pattern_ids))

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

    elapsed_s = (time_ns() - start_time) / 1e9

    # Build SimInfo
    info = SimInfo(
        elapsed_s=elapsed_s,
        backend=:cpu,
        device_id=-1,
        seed=nothing,
        smld_true=smld_true,
        smld_model=smld_model,
        n_patterns=n_patterns,
        n_emitters=length(smld_true.emitters),
        n_localizations=length(smld_noisy.emitters),
        n_frames=params.nframes
    )

    return smld_noisy, info
end

# The largest distance of a pattern's binding sites from its center, the point uniform2D and
# uniform3D place it at and rotate it about (0 for a point pattern). A 3D pattern is rotated in
# 3D, so its full 3D distance bounds how far a site can reach sideways. Internal.
_pattern_radius(p::Pattern2D) = isempty(p.x) ? 0.0 : Float64(maximum(hypot.(p.x, p.y)))
_pattern_radius(p::Pattern3D) = isempty(p.x) ? 0.0 :
    Float64(maximum(sqrt.(p.x .^ 2 .+ p.y .^ 2 .+ p.z .^ 2)))
