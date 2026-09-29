"""
    Analysis functions for diffusion simulations.

This file contains functions for analyzing diffusion simulations,
including dimer detection and analysis of diffusion dynamics.
"""

"""
    get_dimers(smld::BasicSMLD)

Extract a new BasicSMLD containing only emitters in dimer state.

# Arguments
- `smld::BasicSMLD`: Original SMLD with all emitters

# Returns
- `BasicSMLD`: New SMLD containing only dimers

# Example
```julia
# Extract only dimers from simulation results
smld = simulate(params)
dimer_smld = get_dimers(smld)
```
"""
function get_dimers(smld::BasicSMLD)
    # Extract emitters in dimer state
    dimer_emitters = filter(e -> e.state == :dimer, smld.emitters)
    
    # Create new SMLD with same parameters but only dimer emitters
    BasicSMLD(
        dimer_emitters,
        smld.camera,
        smld.n_frames,
        smld.n_datasets,
        copy(smld.metadata)
    )
end

"""
    get_monomers(smld::BasicSMLD)

Extract a new BasicSMLD containing only emitters in monomer state.

# Arguments
- `smld::BasicSMLD`: Original SMLD with all emitters

# Returns
- `BasicSMLD`: New SMLD containing only monomers
"""
function get_monomers(smld::BasicSMLD)
    # Extract emitters in monomer state
    monomer_emitters = filter(e -> e.state == :monomer, smld.emitters)
    
    # Create new SMLD with same parameters but only monomer emitters
    BasicSMLD(
        monomer_emitters,
        smld.camera,
        smld.n_frames,
        smld.n_datasets,
        copy(smld.metadata)
    )
end

"""
    filter_by_state(smld::BasicSMLD, state::Symbol)

Filter emitters by their state (monomer or dimer).

# Arguments
- `smld::BasicSMLD`: The original SMLD with diffusing emitters
- `state::Symbol`: State to filter by (:monomer or :dimer)

# Returns
- `BasicSMLD`: New SMLD containing only emitters with the specified state

# Example
```julia
# Get only monomers
monomer_smld = filter_by_state(smld, :monomer)

# Get only dimers
dimer_smld = filter_by_state(smld, :dimer)
```
"""
function filter_by_state(smld::BasicSMLD, state::Symbol)
    # Filter emitters by state
    filtered_emitters = filter(e -> e.state == state, smld.emitters)
    
    # Create new SMLD with same metadata
    return BasicSMLD(
        filtered_emitters,
        smld.camera,
        smld.n_frames,
        smld.n_datasets,
        copy(smld.metadata)
    )
end

"""
    analyze_dimer_fraction(smld::BasicSMLD)

Calculate the fraction of dimers per frame.

# Arguments
- `smld::BasicSMLD`: SMLD containing all emitters

# Returns
- `Tuple{Vector{Int}, Vector{Float64}}`: Frame numbers and dimer fractions

# Example
```julia
# Calculate dimer fraction over time
frames, fractions = analyze_dimer_fraction(smld)
plot(frames, fractions, xlabel="Frame", ylabel="Dimer Fraction")
```
"""
function analyze_dimer_fraction(smld::BasicSMLD)
    # Get all unique frame numbers
    frames = sort(unique([e.frame for e in smld.emitters]))
    
    # Calculate dimer fraction for each frame
    fractions = Float64[]
    for frame in frames
        # Get emitters in this frame
        frame_emitters = filter(e -> e.frame == frame, smld.emitters)
        
        # Get unique molecule IDs (since dimers have two emitters)
        unique_ids = unique([e.track_id for e in frame_emitters])
        n_molecules = length(unique_ids)
        
        # Count dimers
        dimer_emitters = filter(e -> e.state == :dimer, frame_emitters)
        # Each dimer appears twice, so divide by 2
        n_dimers = length(dimer_emitters) ÷ 2
        
        # Calculate fraction
        fraction = n_molecules > 0 ? n_dimers / n_molecules : 0.0
        push!(fractions, fraction)
    end
    
    return frames, fractions
end

"""
    track_state_changes(smld::BasicSMLD)

Track state changes of molecules over time.

# Arguments
- `smld::BasicSMLD`: SMLD containing all emitters

# Returns
- `Dict{Int, Vector{Tuple{Int, Symbol}}}`: Dictionary mapping molecule IDs to
  vectors of (frame, state) pairs

# Example
```julia
# Track state changes of molecules
state_history = track_state_changes(smld)

# Plot state history for molecule 1
history = state_history[1]
frames = [h[1] for h in history]
states = [h[2] for h in history]
```
"""
function track_state_changes(smld::BasicSMLD)
    # Get all unique molecule IDs
    molecule_ids = unique([e.track_id for e in smld.emitters])
    
    # Initialize state history
    state_history = Dict{Int, Vector{Tuple{Int, Symbol}}}()
    
    # For each molecule, track state changes
    for id in molecule_ids
        # Get all emitters for this molecule, ordered by frame
        mol_emitters = filter(e -> e.track_id == id, smld.emitters)
        sort!(mol_emitters, by = e -> e.frame)
        
        # Extract frame and state
        history = [(e.frame, e.state) for e in mol_emitters]
        
        # Remove consecutive duplicates
        unique_history = Vector{Tuple{Int, Symbol}}()
        last_state = nothing
        
        for (frame, state) in history
            if state != last_state
                push!(unique_history, (frame, state))
                last_state = state
            end
        end
        
        state_history[id] = unique_history
    end
    
    return state_history
end

"""
    analyze_dimer_lifetime(smld::BasicSMLD)

Calculate the average lifetime of dimers.

# Arguments
- `smld::BasicSMLD`: SMLD containing all emitters

# Returns
- `Float64`: Average dimer lifetime in seconds

"""
function analyze_dimer_lifetime(smld::BasicSMLD)
    # Track state changes
    state_history = track_state_changes(smld)
    
    # Calculate dimer lifetimes
    lifetimes = Float64[]
    
    for (id, history) in state_history
        # Find all dimer periods
        dimer_start = nothing
        
        for i in 1:length(history)
            frame, state = history[i]
            
            # Start of dimer period
            if state == :dimer && (i == 1 || history[i-1][2] != :dimer)
                dimer_start = frame
            end
            
            # End of dimer period
            if dimer_start !== nothing && (state != :dimer || i == length(history))
                dimer_end = frame
                
                # Convert frames to time
                params = smld.metadata["simulation_parameters"]
                t_start = (dimer_start - 1) / params.camera_framerate
                t_end = (dimer_end - 1) / params.camera_framerate
                
                # Calculate lifetime
                lifetime = t_end - t_start
                push!(lifetimes, lifetime)
                
                dimer_start = nothing
            end
        end
    end
    
    # Calculate average lifetime
    return isempty(lifetimes) ? 0.0 : mean(lifetimes)
end

"""
    frame_dimer_truth(smld::BasicSMLD)

Per-frame dimer ground truth for every molecule in a diffusion simulation.

Returns a `Vector` of `NamedTuple{(:frame, :track_id, :partner_id, :bound_fraction, :t_form, :t_break, :mixed)}`,
one row per (frame, track_id) present, sorted by (frame, track_id).

- `bound_fraction`: fraction of the track's records in that frame with `state == :dimer`
- `partner_id::Int`: partner of the last bound record in the frame, or `0` if none was bound
- `t_form::Float64`: timestamp of the first record in the frame that is `:dimer` while the
  track's previous record (in time, across frames) was `:monomer`; `NaN` if none
- `t_break::Float64`: the same for `:monomer` after `:dimer`
- `mixed::Bool`: `true` when the row's `partner_id` is nonzero and exactly one of the track and its
  partner is immobile (monomer D == 0), read from `smld.metadata["monomer_D"]`, or for an SMLD without
  a track there from its `"monomer_class"` and the `monomer_mobility` of `"simulation_parameters"`;
  `false` when unbound or when neither D is known. Two mobile partners with different D, or two
  immobile ones, are not mixed.

This is sub-step resolution. With exposure shorter than the frame period, a change that
happens during the gap between exposures shows at the next frame's first record.
"""
function frame_dimer_truth(smld::BasicSMLD)
    by_track = Dict{Int,Vector{Int}}()
    for (i, e) in enumerate(smld.emitters)
        push!(get!(by_track, e.track_id, Int[]), i)
    end

    Dsaved = get(smld.metadata, "monomer_D", nothing)
    class = get(smld.metadata, "monomer_class", nothing)
    cfg = get(smld.metadata, "simulation_parameters", nothing)
    mix = cfg isa DiffusionSMLMConfig ? cfg.monomer_mobility : Tuple{Float64,Float64}[]
    D_of(id) = Dsaved !== nothing && haskey(Dsaved, id) ? Float64(Dsaved[id]) :
               class !== nothing && haskey(class, id) && 1 <= class[id] <= length(mix) ? mix[class[id]][2] : nothing
    function is_mixed(id, partner)
        partner == 0 && return false
        D1, D2 = D_of(id), D_of(partner)
        return D1 !== nothing && D2 !== nothing && ((D1 == 0) ⊻ (D2 == 0))
    end

    rows = NamedTuple{(:frame, :track_id, :partner_id, :bound_fraction, :t_form, :t_break, :mixed),
                      Tuple{Int,Int,Int,Float64,Float64,Float64,Bool}}[]
    for (id, idx) in by_track
        sort!(idx, by = i -> smld.emitters[i].timestamp)
        recs = smld.emitters[idx]
        # Records are time sorted, so each frame is one contiguous run
        n = 0
        n_bound = 0
        partner = 0
        t_form = NaN
        t_break = NaN
        for (k, e) in enumerate(recs)
            if k > 1 && e.frame != recs[k-1].frame
                push!(rows, (frame=recs[k-1].frame, track_id=id, partner_id=partner,
                             bound_fraction=n_bound / n, t_form=t_form, t_break=t_break,
                             mixed=is_mixed(id, partner)))
                n = 0
                n_bound = 0
                partner = 0
                t_form = NaN
                t_break = NaN
            end
            n += 1
            if e.state == :dimer
                n_bound += 1
                partner = something(e.partner_id, 0)
            end
            if k > 1
                prev = recs[k-1].state
                if isnan(t_form) && e.state == :dimer && prev == :monomer
                    t_form = Float64(e.timestamp)
                elseif isnan(t_break) && e.state == :monomer && prev == :dimer
                    t_break = Float64(e.timestamp)
                end
            end
        end
        push!(rows, (frame=recs[end].frame, track_id=id, partner_id=partner,
                     bound_fraction=n_bound / n, t_form=t_form, t_break=t_break,
                     mixed=is_mixed(id, partner)))
    end
    sort!(rows, by = r -> (r.frame, r.track_id))
    return rows
end
