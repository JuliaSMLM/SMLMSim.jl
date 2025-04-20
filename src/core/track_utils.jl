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
    # Filter emitters by track_id directly using Julia's built-in filter function
    # instead of relying on the @filter macro which might have issues
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
    track_ids = Set{Int}()
    for e in smld.emitters
        push!(track_ids, e.track_id)
    end
    
    # Create a new SMLD for each track
    sorted_track_ids = sort(collect(track_ids))
    track_smlds = Vector{BasicSMLD}(undef, length(sorted_track_ids))
    
    for (i, id) in enumerate(sorted_track_ids)
        # Directly filter emitters for this track_id
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