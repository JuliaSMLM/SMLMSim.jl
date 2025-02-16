"""
    CTMC{T<:AbstractFloat, U<:Int}

A Continuous Time Markov Chain representation storing the full trajectory of state transitions.

# Fields
- `simulation_time::T`: Total simulation time span
- `transitiontimes::Vector{T}`: Time points at which state changes occurred, starting at 0.0 
- `states::Vector{U}`: Sequence of states entered at each transition time, starting with initial state

# Type Parameters
- `T`: Floating point type for time values
- `U`: Integer type for state indices

# Note
The `states` and `transitiontimes` vectors have the same length, with each entry in `states[i]`
representing the state entered at time `transitiontimes[i]`. The system remains in `states[i]` 
until time `transitiontimes[i+1]`.
"""
mutable struct CTMC{T<:AbstractFloat,U<:Int}
    simulation_time::T
    transitiontimes::Vector{T}
    states::Vector{U}
end

"""
    CTMC(q::Array{T}, simulation_time::T, state1::Int) where {T<:AbstractFloat}

Construct a Continuous Time Markov Chain simulation from a rate matrix.

# Arguments
- `q::Array{T}`: Rate matrix where q[i,j] is the transition rate from state i to j
- `simulation_time::T`: Total time to simulate
- `state1::Int`: Initial state

# Returns
- `CTMC{T,Int}`: Simulated CTMC with transition times and states

# Details
Simulates a CTMC using the Gillespie algorithm:
1. Start in state1 at time 0
2. For current state i:
   - Calculate total exit rate k_tot = Σ_j q[i,j]
   - Sample time until next transition from Exp(k_tot)
   - Sample next state j with probability q[i,j]/k_tot
3. Repeat until exceeding simulation_time
"""
function CTMC(q::Array{T}, simulation_time::T, state1::Int) where {T<:AbstractFloat}
    lastchange = zero(T)
    currentstate = state1

    states = [state1]
    transitiontimes = [lastchange]

    sidx = Vector((1:size(q, 1)))

    while lastchange < simulation_time
        # get time for state change
        k_tot = sum(q[currentstate, :])
        Δt = rand(Exponential(1 / k_tot))

        # get the new state
        ps = q[currentstate, :]
        ps = ps ./ sum(ps)
        deleteat!(ps, currentstate)
        xs = sidx[:]
        deleteat!(xs, currentstate)
        newstate = rand(DiscreteNonParametric(xs, ps))

        # update CTMC
        lastchange += Δt
        push!(states, newstate)
        push!(transitiontimes, lastchange)
        currentstate = newstate
    end
    return CTMC(simulation_time, transitiontimes, states)
end

"""
    getstate(ctmc::CTMC, t::AbstractFloat)

Get the state of the CTMC at a specific time point.

# Arguments
- `ctmc::CTMC`: The CTMC to query
- `t::AbstractFloat`: Time point of interest

# Returns
- `Int`: State of the chain at time t

# Note
Searches through transition times to find the state active at time t.
Returns the state that was entered at the last transition before t.
"""
function getstate(ctmc::CTMC, t::AbstractFloat)
    nstates = length(ctmc.states)
    for nn = 2:nstates
        if t < ctmc.transitiontimes[nn]
            return ctmc.states[nn-1]
        end
    end
end

"""
    getnext(ctmc::CTMC, t::AbstractFloat)

Get the next state transition after a specific time point.

# Arguments
- `ctmc::CTMC`: The CTMC to query
- `t::AbstractFloat`: Current time point

# Returns
- `Tuple{Int,AbstractFloat}`: (next_state, transition_time)

# Note
Returns the next state that will be entered and when it will be entered,
searching from the current time point forward.
"""
function getnext(ctmc::CTMC, t::AbstractFloat)
    nstates = length(ctmc.states)
    for nn = 1:nstates
        if t < ctmc.transitiontimes[nn]
            return ctmc.states[nn], ctmc.transitiontimes[nn]
        end
    end
end

"""
    intensitytrace(f::GenericFluor, nframes::Int, framerate::Real; state1=1)

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

# Note
- State 1 is assumed to be the fluorescent state
- Emission only occurs in state 1 with rate f.γ
- Frame exposure is assumed to be 1/framerate (100% duty cycle)
"""
function intensitytrace(f::GenericFluor, nframes::Int, framerate::Real; state1=1)
    simulation_time = nframes / framerate

    # generate CTMC
    ctmc = SMLMSim.CTMC(f.q, simulation_time, state1)

    # generate integrated photons 
    exptime = 1 / framerate
    photons = zeros(nframes)
    for nn = 1:nframes
        t = (nn - 1) * exptime
        frameend = nn * exptime
        while t < frameend
            currentstate = SMLMSim.getstate(ctmc, t)
            (_, nexttime) = SMLMSim.getnext(ctmc, t)
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
    kineticmodel(smld::BasicSMLD, f::Molecule, nframes::Int, framerate::Real;
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

# Note
The emitter type (2D/3D) is automatically determined from the input SMLD.
Position uncertainties are initialized to 0 and can be set using the
noise() function.
"""
function kineticmodel(smld::BasicSMLD, f::Molecule, nframes::Int, framerate::Real;
                     ndatasets::Int=1, minphotons=50.0, state1::Int=2)
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
        photons = SMLMSim.intensitytrace(f, nframes, framerate; state1=state1)
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

"""
    noise(smld::BasicSMLD, σ_psf::Union{AbstractFloat, Vector{<:AbstractFloat}})

Add localization uncertainty to emitter positions based on photon counts.

# Arguments
- `smld::BasicSMLD`: Input SMLD containing emitters
- `σ_psf::Union{AbstractFloat, Vector{<:AbstractFloat}}`: PSF width(s) in microns
    - For 2D: single σ value
    - For 3D: vector [σx, σy, σz]

# Returns
- `BasicSMLD`: New SMLD with noisy positions and updated uncertainties

# Details
For each emitter:
1. Calculates position uncertainty as σ_psf/√N where N is photon count
2. Adds Gaussian noise to positions with appropriate σ
3. Updates uncertainty fields in emitter
4. Preserves all other emitter properties (frame, dataset, track_id)

# Note
Automatically handles both 2D and 3D cases based on emitter type in SMLD.
Input σ_psf must match dimensionality (scalar for 2D, vector for 3D).
"""
function noise(smld::BasicSMLD, σ_psf::Union{AbstractFloat, Vector{<:AbstractFloat}})
    etype = eltype(smld.emitters)
    is_3d = etype <: Emitter3DFit
    
    if is_3d && !isa(σ_psf, Vector)
        error("3D emitter type requires vector of σ_psf values")
    elseif !is_3d && !isa(σ_psf, AbstractFloat)
        error("2D emitter type requires scalar σ_psf value")
    end
    
    new_emitters = similar(smld.emitters)
    
    for (i, emitter) in enumerate(smld.emitters)
        σ = is_3d ? σ_psf ./ sqrt(emitter.photons) : σ_psf / sqrt(emitter.photons)
        
        coords = if is_3d
            (
                emitter.x + randn() * σ[1],
                emitter.y + randn() * σ[2],
                emitter.z + randn() * σ[3]
            )
        else
            (
                emitter.x + randn() * σ,
                emitter.y + randn() * σ
            )
        end
        
        uncertainties = if is_3d
            (σ[1], σ[2], σ[3])
        else
            (σ, σ)
        end
        
        common_params = (
            emitter.photons,
            emitter.bg,
            emitter.σ_photons,
            emitter.σ_bg
        )
        
        new_emitters[i] = etype(
            coords...,            # Position coordinates
            common_params...,     # Photons and background
            uncertainties...,     # Position uncertainties
            frame=emitter.frame,
            dataset=emitter.dataset,
            track_id=emitter.track_id,
            id=emitter.id
        )
    end
    
    return BasicSMLD(
        new_emitters,
        smld.camera,
        smld.n_frames,
        smld.n_datasets,
        copy(smld.metadata)
    )
end