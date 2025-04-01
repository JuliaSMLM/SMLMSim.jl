"""
    Photophysics Module

This module provides functions for modeling the photophysics of fluorescent molecules,
including blinking kinetics, intensity traces, and related phenomena.
"""


"""
    intensity_trace(f::GenericFluor, nframes::Int, framerate::Real; state1=1)

Calculate a fluorescence intensity trace by integrating emission during fluorescent state occupancy.

# Arguments
- `f::GenericFluor`: Fluorophore model containing transition rates (q) and emission rate (γ)
- `nframes::Int`: Number of frames to simulate
- `framerate::Real`: Frame rate in Hz
- `state1::Int=1`: Initial state (default: 1 for fluorescent state)

# Returns
- `Vector{Float64}`: Integrated photon counts for each frame

# Details
For each frame:
1. Determines state occupancy using CTMC
2. Integrates emission (rate f.γ) during fluorescent state periods
3. Accumulates photons within frame exposure time (1/framerate)

# Example
```julia
fluor = GenericFluor(; γ=10000.0, q=[0 10; 1e-1 0])
photons = intensity_trace(fluor, 1000, 10.0)
```

# Note
- State 1 is assumed to be the fluorescent state
- Emission only occurs in state 1 with rate f.γ
- Frame exposure is assumed to be 1/framerate (100% duty cycle)
"""
function intensity_trace(f::GenericFluor, nframes::Int, framerate::Real; state1=1)
    # Input validation
    if nframes <= 0
        throw(ArgumentError("Number of frames must be positive"))
    end
    if framerate <= 0
        throw(ArgumentError("Frame rate must be positive"))
    end
    
    simulation_time = nframes / framerate

    # generate CTMC
    ctmc = CTMC(f.q, simulation_time, state1)

    # generate integrated photons 
    exptime = 1 / framerate
    photons = zeros(nframes)
    for nn = 1:nframes
        t = (nn - 1) * exptime
        frameend = nn * exptime
        while t < frameend
            currentstate = get_state(ctmc, t)
            (_, nexttime) = get_next(ctmc, t)
            tend = min(frameend, nexttime)
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
                 ndatasets::Int=1, minphotons=50.0, state1::Int=2)

Generate kinetic blinking model from existing localization data.

# Arguments
- `smld::BasicSMLD`: Input SMLD containing true emitter positions
- `f::Molecule`: Fluorophore model with kinetic rates
- `nframes::Int`: Number of frames to simulate
- `framerate::Real`: Frame rate in Hz
- `ndatasets::Int=1`: Number of independent datasets to generate
- `minphotons::Float64=50.0`: Minimum photons for detection
- `state1::Int=2`: Initial state (default: 2 for dark state)

# Returns
- `BasicSMLD`: New SMLD with simulated blinking kinetics

# Details
For each unique position in the input SMLD:
1. Simulates fluorophore blinking using the kinetic model
2. Creates emitters for frames where photon count exceeds threshold
3. Preserves track_id for linking emitters from same position
4. Maintains camera and extends metadata from input SMLD

# Example
```julia
camera = IdealCamera(1:128, 1:128, 0.1)
pattern = Nmer2D()
smld_true, _, _ = simulate(pattern=pattern, camera=camera)

# Add blinking kinetics
fluor = GenericFluor(; γ=10000.0, q=[0 10; 1e-1 0])
smld_model = kinetic_model(smld_true, fluor, 1000, 10.0)
```

# Note
The emitter type (2D/3D) is automatically determined from the input SMLD.
Position uncertainties are initialized to 0 and can be set using the
apply_noise() function.
"""
function kinetic_model(smld::BasicSMLD, f::Molecule, nframes::Int, framerate::Real;
                      ndatasets::Int=1, minphotons=50.0, state1::Int=2)
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
    
    emitter_type = eltype(smld.emitters)
    emitters = Vector{emitter_type}()

    # Get true positions from input SMLD
    grouped_emitters = Dict()
    for e in smld.emitters
        key = e.track_id
        if !haskey(grouped_emitters, key)
            grouped_emitters[key] = e
        end
    end

    # Sort by track_id to maintain order
    true_positions = sort(collect(values(grouped_emitters)), by=e -> e.track_id)

    for dd = 1:ndatasets, pos in true_positions
        photons = intensity_trace(f, nframes, framerate; state1=state1)
        framenum = findall(photons .> minphotons)

        # Create emitter for each frame where photons > threshold
        for frame in framenum
            if emitter_type <: Emitter2DFit
                emitter = Emitter2DFit{Float64}(
                    pos.x, pos.y,        # positions
                    photons[frame],      # photons
                    0.0,                 # background
                    0.0, 0.0,            # σ_x, σ_y
                    0.0, 0.0;            # σ_photons, σ_bg
                    frame=frame,
                    dataset=dd,
                    track_id=pos.track_id
                )
            elseif emitter_type <: Emitter3DFit
                emitter = Emitter3DFit{Float64}(
                    pos.x, pos.y, pos.z,  # positions
                    photons[frame],       # photons
                    0.0,                  # background
                    0.0, 0.0, 0.0,        # σ_x, σ_y, σ_z
                    0.0, 0.0;             # σ_photons, σ_bg
                    frame=frame,
                    dataset=dd,
                    track_id=pos.track_id
                )
            else
                error("Unsupported emitter type: $emitter_type")
            end
            push!(emitters, emitter)
        end
    end

    # Create new metadata combining original and simulation parameters
    metadata = copy(smld.metadata)
    metadata["simulation_type"] = "kinetic_model"
    metadata["framerate"] = framerate

    return BasicSMLD(emitters, smld.camera, nframes, ndatasets, metadata)
end
