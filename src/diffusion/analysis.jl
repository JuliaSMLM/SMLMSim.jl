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
        unique_ids = unique([e.id for e in frame_emitters])
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
    molecule_ids = unique([e.id for e in smld.emitters])
    
    # Initialize state history
    state_history = Dict{Int, Vector{Tuple{Int, Symbol}}}()
    
    # For each molecule, track state changes
    for id in molecule_ids
        # Get all emitters for this molecule, ordered by frame
        mol_emitters = filter(e -> e.id == id, smld.emitters)
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
