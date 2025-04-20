"""
    Track utility functions for SMLD objects.

This file contains helper functions for working with tracks (trajectories) in SMLD data.
"""

"""
    get_track(smld::BasicSMLD, id::Int)

Return a new SMLD containing only emitters with the specified track_id.

# Arguments
- `smld::BasicSMLD`: The original SMLD
- `id::Int`: Track ID to filter by

# Returns
- `BasicSMLD`: New SMLD containing only emitters from the specified track

# Example
```julia
# Get all emitters belonging to track 5
track_smld = get_track(smld, 5)
```
"""
function get_track(smld::BasicSMLD, id::Int)
    # Filter emitters by track_id
    track_emitters = @filter(smld, track_id == id)
    
    # Create new SMLD with same metadata
    return BasicSMLD(
        track_emitters.emitters,
        smld.camera,
        smld.n_frames,
        smld.n_datasets,
        copy(smld.metadata)
    )
end

"""
    get_num_tracks(smld::BasicSMLD)

Return the number of unique tracks (based on track_id) in the SMLD.

# Arguments
- `smld::BasicSMLD`: The SMLD to analyze

# Returns
- `Int`: Number of unique track IDs

# Example
```julia
# Get the number of tracks
n_tracks = get_num_tracks(smld)
```
"""
function get_num_tracks(smld::BasicSMLD)
    # Use a Set for efficient uniqueness checking
    track_ids = Set{Int}()
    
    for e in smld.emitters
        push!(track_ids, e.track_id)
    end
    
    return length(track_ids)
end

"""
    get_tracks(smld::BasicSMLD)

Return a vector of SMLD objects, one for each unique track (based on track_id).

# Arguments
- `smld::BasicSMLD`: The original SMLD

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
function get_tracks(smld::BasicSMLD)
    # Group emitters by track_id more efficiently
    emitters_by_track = Dict{Int, Vector{eltype(smld.emitters)}}()
    
    for e in smld.emitters
        if !haskey(emitters_by_track, e.track_id)
            emitters_by_track[e.track_id] = eltype(smld.emitters)[]
        end
        push!(emitters_by_track[e.track_id], e)
    end
    
    # Create a new SMLD for each track
    track_ids = sort(collect(keys(emitters_by_track)))
    track_smlds = Vector{BasicSMLD}(undef, length(track_ids))
    
    for (i, id) in enumerate(track_ids)
        track_smlds[i] = BasicSMLD(
            emitters_by_track[id],
            smld.camera,
            smld.n_frames,
            smld.n_datasets,
            copy(smld.metadata)
        )
    end
    
    return track_smlds
end