"""
    Track utility functions for diffusing emitters in SMLD objects.

Specialized version of track utility functions for diffusing emitters.
We now use track_id consistently across all emitter types, so these functions
are maintained for backward compatibility but will eventually be deprecated.
"""

"""
    get_track(smld::BasicSMLD, id::Int)

Return a new SMLD containing only emitters with the specified track_id.
Specialized version for diffusing emitters, maintained for backward compatibility.

# Arguments
- `smld::BasicSMLD`: The original SMLD with diffusing emitters
- `id::Int`: track_id to filter by

# Returns
- `BasicSMLD`: New SMLD containing only emitters from the specified track

# Example
```julia
# Get all emitters belonging to ID 5
track_smld = get_track(smld, 5)
```
"""
function get_track(smld::BasicSMLD{T, E}, id::Int) where {T, E <: AbstractDiffusingEmitter}
    # Filter emitters by track_id directly
    filtered_emitters = filter(e -> e.track_id == id, smld.emitters)
    
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
    get_num_tracks(smld::BasicSMLD)

Return the number of unique tracks (based on track_id) in the SMLD with diffusing emitters.
Specialized version for diffusing emitters, maintained for backward compatibility.

# Arguments
- `smld::BasicSMLD`: The SMLD with diffusing emitters to analyze

# Returns
- `Int`: Number of unique IDs

# Example
```julia
# Get the number of tracks
n_tracks = get_num_tracks(smld)
```
"""
function get_num_tracks(smld::BasicSMLD{T, E}) where {T, E <: AbstractDiffusingEmitter}
    # Use a Set for efficient uniqueness checking
    unique_ids = Set{Int}()
    
    for e in smld.emitters
        push!(unique_ids, e.track_id)
    end
    
    return length(unique_ids)
end

"""
    get_tracks(smld::BasicSMLD)

Return a vector of SMLD objects, one for each unique track (based on track_id) for diffusing emitters.
Specialized version for diffusing emitters, maintained for backward compatibility.

# Arguments
- `smld::BasicSMLD`: The original SMLD with diffusing emitters

# Returns
- `Vector{BasicSMLD}`: Vector of SMLD objects, one per track

# Example
```julia
# Get all tracks as separate SMLD objects
track_smlds = get_tracks(smld)

# Access the first track
first_track = track_smlds[1]
```
"""
function get_tracks(smld::BasicSMLD{T, E}) where {T, E <: AbstractDiffusingEmitter}
    # Group emitters by track_id
    unique_ids = Set{Int}()
    for e in smld.emitters
        push!(unique_ids, e.track_id)
    end
    
    # Create a new SMLD for each track
    sorted_ids = sort(collect(unique_ids))
    track_smlds = Vector{BasicSMLD}(undef, length(sorted_ids))
    
    for (i, id) in enumerate(sorted_ids)
        # Directly filter emitters for this id
        filtered_emitters = filter(e -> e.track_id == id, smld.emitters)
        
        track_smlds[i] = BasicSMLD(
            filtered_emitters,
            smld.camera,
            smld.n_frames,
            smld.n_datasets,
            copy(smld.metadata)
        )
    end
    
    return track_smlds
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
function filter_by_state(smld::BasicSMLD{T, E}, state::Symbol) where {T, E <: AbstractDiffusingEmitter}
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
