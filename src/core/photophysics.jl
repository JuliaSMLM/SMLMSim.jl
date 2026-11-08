"""
    Photophysics Module

This module provides functions for modeling the photophysics of fluorescent molecules,
including blinking kinetics, intensity traces, and sampling from equilibrium distributions.

# Make sure LinearAlgebra and Distributions are imported at the Core module level
# for norm() and rand() functions
"""


"""
    intensity_trace(f::GenericFluor, nframes::Int, framerate::Real; state1=1, burn_in=0.0)

Calculate a fluorescence intensity trace by integrating emission during fluorescent state occupancy.

# Arguments
- `f::GenericFluor`: Fluorophore model containing transition rates (q) and emission rate (γ)
- `nframes::Int`: Number of frames to simulate
- `framerate::Real`: Frame rate in Hz
- `state1::Int=1`: Initial state (default: 1 for fluorescent state)
- `burn_in::Real=0.0`: Pre-illumination time in seconds before recording starts.
  The CTMC runs for this duration first, allowing the system to reach a pseudo-equilibrium
  (e.g., some molecules bleach before data collection begins).

# Returns
- `Vector{Float64}`: Integrated photon counts for each frame

# Details
For each frame:
1. Determines state occupancy using CTMC (which starts at time 0, but recording starts at burn_in)
2. Integrates emission (rate f.γ) during fluorescent state periods
3. Accumulates photons within frame exposure time (1/framerate)

# Example
```julia
fluor = GenericFluor(10000.0, [-10.0 10.0; 1e-1 -1e-1])
photons = intensity_trace(fluor, 1000, 10.0)

# With 5 second burn-in (pre-illumination before recording)
photons = intensity_trace(fluor, 1000, 10.0; burn_in=5.0)
```

# Note
- State 1 is assumed to be the fluorescent state
- Emission only occurs in state 1 with rate f.γ
- Frame exposure is assumed to be 1/framerate (100% duty cycle)
- burn_in is useful for models with photobleaching to simulate pre-illumination
- Through 0.7.3 a bleached molecule could hang when the final frame end rounded past the
  CTMC's end time. A state with no transition left now holds through the frame end.
- The rate matrix and simulation time are promoted to a common floating point type when
  building the chain. Through 0.7.3 a Float32 rate matrix with a Float64 time failed with a
  method error. Float64 simulations keep the same random draws and photon counts.
"""
function intensity_trace(f::GenericFluor, nframes::Int, framerate::Real; state1=1, burn_in::Real=0.0)
    # Input validation
    if nframes <= 0
        throw(ArgumentError("Number of frames must be positive"))
    end
    if framerate <= 0
        throw(ArgumentError("Frame rate must be positive"))
    end
    if burn_in < 0
        throw(ArgumentError("burn_in must be non-negative"))
    end

    # Total simulation time includes burn-in period
    recording_time = nframes / framerate
    simulation_time = burn_in + recording_time

    # generate CTMC for full duration (burn_in + recording)
    T = promote_type(eltype(f.q), typeof(simulation_time))
    q = eltype(f.q) === T ? f.q : T.(f.q)
    ctmc = CTMC(q, T(simulation_time), state1)

    # generate integrated photons (recording starts after burn_in)
    exptime = 1 / framerate
    photons = zeros(nframes)
    for nn = 1:nframes
        t = burn_in + (nn - 1) * exptime  # offset by burn_in
        frameend = burn_in + nn * exptime
        while t < frameend
            currentstate = get_state(ctmc, T(t))
            (_, nexttime) = get_next(ctmc, T(t))
            tend = min(frameend, nexttime)
            tend <= t && (tend = frameend)   # no transition left: the state holds to the frame end
            if currentstate == 1  # the fluorescent state
                photons[nn] += f.γ * (tend - t)
            end
            t = tend
        end
    end
    return photons
end

"""
    kinetic_model(smld::BasicSMLD, f::Molecule, nframes::Int, framerate::Real;
                 ndatasets::Int=1, minphotons=50.0, state1::Union{Int, Symbol}=:equilibrium,
                 burn_in::Real=0.0,
                 brightness::Union{Nothing, Distributions.UnivariateDistribution}=nothing)

Generate kinetic blinking model from existing localization data.

# Arguments
- `smld::BasicSMLD`: Input SMLD containing true emitter positions
- `f::Molecule`: Fluorophore model with kinetic rates
- `nframes::Int`: Number of frames to simulate
- `framerate::Real`: Frame rate in Hz
- `ndatasets::Int=1`: Number of independent datasets to generate
- `minphotons::Float64=50.0`: Minimum photons for detection
- `state1::Union{Int, Symbol}=:equilibrium`: Initial state specification:
  - `::Int`: Specific state to start in (1=on, 2=off typically)
  - `:equilibrium`: Sample from equilibrium distribution (default)
- `burn_in::Real=0.0`: Pre-illumination time in seconds before recording starts.
  Simulates the common experimental protocol where high laser power is applied
  for several seconds before data collection begins, allowing some molecules
  to bleach and others to reach a pseudo-equilibrium blinking state.
- `brightness::Union{Nothing, Distributions.UnivariateDistribution}=nothing`: Distribution of a
  dimensionless factor `b` on the molecule model's emission rate. A molecule with factor `b`
  emits at `b * f.γ` photons/s for all of its frames; the rate matrix is the model's, unchanged.
  `b` is drawn once per molecule and dataset, just before that molecule's initial state, and
  converted to `Float64`. Redrawing per dataset is intended: each dataset is an independent
  sample, as its initial states and blinking already are, so the same position can differ in
  brightness between datasets. The factor must be finite and not negative. Before any random
  draw, `ArgumentError` is thrown when `minimum(brightness)` is negative or NaN, or when
  `brightness` is given with a model that is not a `GenericFluor`. Every drawn factor is also
  checked: a draw that is NaN, infinite or negative throws `ArgumentError` naming the value and
  the distribution (a support can hold such a value above a valid minimum, as in
  `DiscreteNonParametric([1.0, NaN], [0.5, 0.5])`). To keep a normal distribution at or above
  a floor, use either `truncated(d; lower=floor)` (cut at the floor and renormalized; no draw
  lies at the floor) or `censored(d; lower=floor)` (draws below the floor are set to it).
  `nothing` (default) gives every molecule the rate `f.γ`, with no extra random draw.

# Returns
- `BasicSMLD`: New SMLD with simulated blinking kinetics. When `brightness` is given, its metadata
  holds `"brightness"` (the distribution) and `"brightness_factors"`, a
  `Dict{Tuple{Int,Int},Float64}` from `(dataset, track_id)` to the drawn factor, with the
  output `track_id`. Without `brightness` neither key is present: none is added, and any carried
  in from the input's metadata is removed.

# Details
Each non-zero `track_id` identifies one molecule, using its first record. Repeated non-zero
ids at different positions warn once because later positions are dropped. Each emitter with
`track_id == 0` is an independent, unlinked molecule; its output id is the next free integer
above the largest non-zero id, in input order. If these ids would overflow `Int`, an
`ArgumentError` names the largest non-zero id. Molecules run in sorted non-zero id order,
followed by unlinked emitters in input order. Through 0.7.3 all unlinked emitters were grouped
as one molecule, and conflicting positions for a linked id were dropped silently.

For each molecule in the input SMLD:
1. Simulates fluorophore blinking using the kinetic model
2. Creates emitters for frames where photon count exceeds threshold
3. Preserves non-zero track_id and assigns ids to unlinked emitters
4. Maintains camera and extends metadata from input SMLD

# Example
```julia
camera = IdealCamera(1:128, 1:128, 0.1)
pattern = Nmer2D()
_, info = simulate(StaticSMLMConfig(); pattern=pattern, camera=camera)
smld_true = info.smld_true

# Add blinking kinetics
fluor = GenericFluor(10000.0, [-10.0 10.0; 1e-1 -1e-1])
smld_model = kinetic_model(smld_true, fluor, 1000, 10.0)

# With 5 second burn-in to simulate pre-illumination
smld_model = kinetic_model(smld_true, fluor, 1000, 10.0; burn_in=5.0)

# Molecule-to-molecule brightness: mean 1, sd 0.3, no factor below 0.7
# (Normal, truncated and censored come from Distributions.jl: `using Distributions`)
smld_model = kinetic_model(smld_true, fluor, 1000, 10.0;
                           brightness=truncated(Normal(1, 0.3); lower=0.7))
smld_model = kinetic_model(smld_true, fluor, 1000, 10.0;
                           brightness=censored(Normal(1, 0.3); lower=0.7))
```

# Note
The emitter type (2D/3D fit) is determined from the individual input emitters. Abstract
container element types and mixed numeric precisions within one dimension are supported;
output fits use Float64 coordinates. Mixed 2D and 3D inputs throw `ArgumentError`.
Position uncertainties are initialized to 0 and can be set using the
apply_noise() function. For models with photobleaching (absorbing states),
use `state1=1` instead of `:equilibrium` since absorbing states have no
true equilibrium.
"""
function kinetic_model(smld::BasicSMLD, f::Molecule, nframes::Int, framerate::Real;
                      ndatasets::Int=1, minphotons=50.0, state1::Union{Int, Symbol}=:equilibrium,
                      burn_in::Real=0.0,
                      brightness::Union{Nothing, Distributions.UnivariateDistribution}=nothing)
    # Input validation
    if nframes <= 0
        throw(ArgumentError("Number of frames must be positive"))
    end
    if framerate <= 0
        throw(ArgumentError("Frame rate must be positive"))
    end
    if ndatasets <= 0
        throw(ArgumentError("Number of datasets must be positive"))
    end
    if minphotons < 0
        throw(ArgumentError("Minimum photons must be non-negative"))
    end
    if burn_in < 0
        throw(ArgumentError("burn_in must be non-negative"))
    end
    _check_brightness(brightness, f)

    # Choose concrete Float64 output storage from the input values, not the container type.
    emitters = if isempty(smld.emitters)
        similar(smld.emitters, 0)
    elseif all(e -> e isa Emitter2DFit, smld.emitters)
        Emitter2DFit{Float64}[]
    elseif all(e -> e isa Emitter3DFit, smld.emitters)
        Emitter3DFit{Float64}[]
    else
        throw(ArgumentError("kinetic_model requires only Emitter2DFit or only Emitter3DFit emitters; mixed 2D and 3D inputs are not supported"))
    end

    # Get true positions from input SMLD
    grouped_emitters = Dict()
    unlinked_emitters = AbstractEmitter[]
    for e in smld.emitters
        key = e.track_id
        if key == 0
            push!(unlinked_emitters, e)
        elseif !haskey(grouped_emitters, key)
            grouped_emitters[key] = e
        else
            first_emitter = grouped_emitters[key]
            if e.x != first_emitter.x || e.y != first_emitter.y ||
               (hasproperty(e, :z) && e.z != first_emitter.z)
                @warn "kinetic_model: repeated non-zero track_id $key has records at different positions; using its first record and dropping later positions" maxlog=1
            end
        end
    end

    # Keep the linked molecule order and append independent, unlinked molecules.
    true_positions = sort(collect(values(grouped_emitters)), by=e -> e.track_id)
    n_linked = length(true_positions)
    largest_id = maximum(keys(grouped_emitters); init=0)
    # Reject exhausted ids before sampling any molecule's initial state or blinking trace.
    try
        Base.Checked.checked_add(largest_id, length(unlinked_emitters))
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("Unlinked track_ids would overflow Int above the largest non-zero id $largest_id"))
    end
    append!(true_positions, unlinked_emitters)

    # Precompute equilibrium distribution if needed (do once, not per emitter)
    initial_state_dist = if state1 == :equilibrium
        q = f.q
        if size(q, 1) == 2  # Common 2-state case (on/off)
            # For 2-state case with negative diagonals
            k_on = q[2, 1]      # off->on rate
            k_off = q[1, 2]     # on->off rate
            # Verify proper structure of rate matrix
            if q[1,1] > 0 || q[2,2] > 0
                @warn "Diagonal elements should be negative; equilibrium calculation may be incorrect"
            end
            
            # Equilibrium probabilities
            p_on = k_on / (k_on + k_off)  # Probability of being in ON state
            [p_on, 1.0 - p_on]  # Distribution vector [p_state1, p_state2]
        else
            # For n-state case, use general algorithm
            compute_equilibrium_distribution(q)
        end
    else
        nothing  # Not needed for explicit state
    end

    brightness_factors = Dict{Tuple{Int,Int},Float64}()
    for dd = 1:ndatasets, (i, pos) in enumerate(true_positions)
        track_id = pos.track_id == 0 ? Base.Checked.checked_add(largest_id, i - n_linked) : pos.track_id
        # One brightness factor per molecule and dataset, drawn before the initial state
        fm = f
        if brightness !== nothing
            b = _draw_brightness(brightness)
            brightness_factors[(dd, track_id)] = b
            fm = GenericFluor(b * f.γ, f.q)
        end
        # Determine initial state based on the option provided
        initial_state = if state1 == :equilibrium
            # Use pre-computed distribution for sampling
            sample_discrete(initial_state_dist)
        else
            # Use specified initial state
            state1
        end
        
        photons = intensity_trace(fm, nframes, framerate; state1=initial_state, burn_in=burn_in)
        framenum = findall(photons .> minphotons)

        # Create emitter for each frame where photons > threshold
        for frame in framenum
            if pos isa Emitter2DFit
                emitter = Emitter2DFit{Float64}(
                    Float64(pos.x), Float64(pos.y),  # positions
                    photons[frame],      # photons
                    0.0,                 # background
                    0.0, 0.0,            # σ_x, σ_y
                    0.0, 0.0;            # σ_photons, σ_bg
                    σ_xy=0.0,            # x-y covariance (0 for symmetric PSF)
                    frame=frame,
                    dataset=dd,
                    track_id=track_id,
                    id=pos.id
                )
            elseif pos isa Emitter3DFit
                emitter = Emitter3DFit{Float64}(
                    Float64(pos.x), Float64(pos.y), Float64(pos.z),  # positions
                    photons[frame],       # photons
                    0.0,                  # background
                    0.0, 0.0, 0.0,        # σ_x, σ_y, σ_z
                    0.0, 0.0;             # σ_photons, σ_bg
                    σ_xy=0.0,             # x-y covariance (0 for symmetric PSF)
                    frame=frame,
                    dataset=dd,
                    track_id=track_id,
                    id=pos.id
                )
            else
                error("Unsupported emitter type: $(typeof(pos))")
            end
            push!(emitters, emitter)
        end
    end

    # Create new metadata combining original and simulation parameters
    metadata = copy(smld.metadata)
    metadata["simulation_type"] = "kinetic_model"
    metadata["framerate"] = framerate
    metadata["initial_state"] = string(state1)
    metadata["burn_in"] = burn_in
    if brightness !== nothing
        metadata["brightness"] = brightness
        metadata["brightness_factors"] = brightness_factors
    else
        # An input from an earlier brightness run must not carry its factors into this one
        delete!(metadata, "brightness")
        delete!(metadata, "brightness_factors")
    end

    return BasicSMLD(emitters, smld.camera, nframes, ndatasets, metadata)
end

# Checks for kinetic_model's and static simulate's `brightness` keyword. Internal.
function _check_brightness(brightness, f::Molecule)
    brightness === nothing && return nothing
    f isa GenericFluor ||
        throw(ArgumentError("brightness requires a GenericFluor molecule model, got $(typeof(f))"))
    # `>=` also rejects a NaN minimum
    minimum(brightness) >= 0 ||
        throw(ArgumentError("brightness must not take negative values (a factor below zero would make a negative emission rate); use truncated(d; lower=floor) to cut the distribution at a floor or censored(d; lower=floor) to clamp draws to it"))
    return nothing
end

# One brightness factor, drawn and checked. The check above sees only the distribution's
# minimum, so a support that also holds NaN or Inf (DiscreteNonParametric([1.0, NaN], ...))
# passes it. One random draw, as an unchecked draw makes. Internal.
function _draw_brightness(brightness)
    x = rand(brightness)
    b = Float64(x)
    isfinite(b) && b >= 0 ||
        throw(ArgumentError("brightness drew the factor $(repr(x)) from $brightness; every drawn factor must be finite and non-negative as a Float64"))
    return b
end

"""
    compute_equilibrium_distribution(q::Matrix{<:AbstractFloat})

Calculate the equilibrium probability distribution for a CTMC rate matrix.

# Arguments
- `q::Matrix{<:AbstractFloat}`: Rate matrix where q[i,j] for i≠j is the transition rate from state i to j,
  and q[i,i] is the negative exit rate from state i

# Returns
- `Vector{Float64}`: Equilibrium probabilities for each state

# Details
For a rate matrix Q, the equilibrium distribution π satisfies π·Q = 0 subject to Σπ = 1.
This function solves the linear system directly to find the equilibrium distribution.
"""
function compute_equilibrium_distribution(q::Matrix{<:AbstractFloat})
    n_states = size(q, 1)
    
    # Check for proper rate matrix structure (negative diagonals or zero for absorbing, rows sum to zero)
    row_sums = sum(q, dims=2)
    diag_vals = diag(q)
    if any(diag_vals .> 0)
        @warn "Some diagonal elements of rate matrix are positive (should be negative or zero)"
    end
    if !all(abs.(row_sums) .< 1e-10)
        @warn "Rate matrix rows do not sum to zero (maximum deviation: $(maximum(abs.(row_sums))))"
    end
    
    # Create the linear system A⋅π = b for equilibrium distribution
    # Replace last row with constraint that probabilities sum to 1
    A = copy(q')
    A[end, :] .= 1.0
    
    # Right-hand side: all zeros except last element = 1
    b = zeros(n_states)
    b[end] = 1.0
    
    # Solve the system
    π = A \ b
    
    # Ensure probabilities are valid (non-negative and sum to 1)
    # This handles numerical precision issues
    π = max.(π, 0)  # Ensure non-negative
    π = π ./ sum(π)  # Normalize to sum to 1
    
    return π
end

"""
    sample_discrete(p::Vector{<:AbstractFloat})

Sample from a discrete probability distribution.

# Arguments
- `p::Vector{<:AbstractFloat}`: Probability distribution

# Returns
- `Int`: Sampled state index

# Details
Samples a state index i with probability p[i] using the inverse CDF method.
"""
function sample_discrete(p::Vector{<:AbstractFloat})
    u = rand()
    cdf = 0.0
    for i in 1:length(p)
        cdf += p[i]
        if u <= cdf
            return i
        end
    end
    return length(p)  # Fallback for numerical precision issues
end
