"""
    Coordinate Noise Module

This module provides functions for modeling the localization coordinate noise
in SMLM, including detection uncertainty calculations based on photon statistics.
"""

"""
    apply_noise(smld::BasicSMLD, σ_psf::Union{AbstractFloat, Vector{<:AbstractFloat}})

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

# Example
```julia
# First create ground truth model with blinking
camera = IdealCamera(1:128, 1:128, 0.1)
pattern = Nmer2D()
_, smld_model, _ = simulate(pattern=pattern, camera=camera)

# Then add localization noise with specific PSF width
smld_noisy = apply_noise(smld_model, 0.13)  # 130nm PSF width
```

# Note
Automatically handles both 2D and 3D cases based on emitter type in SMLD.
Input σ_psf must match dimensionality (scalar for 2D, vector for 3D).
"""
function apply_noise(smld::BasicSMLD, σ_psf::Union{AbstractFloat, Vector{<:AbstractFloat}})
    etype = eltype(smld.emitters)
    is_3d = etype <: Emitter3DFit
    
    if is_3d && !isa(σ_psf, Vector)
        error("3D emitter type requires vector of σ_psf values")
    elseif !is_3d && !isa(σ_psf, AbstractFloat)
        error("2D emitter type requires scalar σ_psf value")
    end
    
    if is_3d && length(σ_psf) != 3
        error("3D emitter type requires vector of 3 σ_psf values [σx, σy, σz]")
    end
    
    new_emitters = similar(smld.emitters)
    
    for (i, emitter) in enumerate(smld.emitters)
        if emitter.photons <= 0
            error("Emitter photon count must be positive")
        end
        
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
        
        # Create new emitter with correct parameter order
        if is_3d
            new_emitters[i] = etype(
                coords[1], coords[2], coords[3],  # x, y, z
                emitter.photons,                  # photons
                emitter.bg,                       # background
                σ[1], σ[2], σ[3],                 # σ_x, σ_y, σ_z
                emitter.σ_photons, emitter.σ_bg;  # σ_photons, σ_bg
                frame=emitter.frame,
                dataset=emitter.dataset,
                track_id=emitter.track_id,
                id=emitter.id
            )
        else
            new_emitters[i] = etype(
                coords[1], coords[2],             # x, y
                emitter.photons,                  # photons
                emitter.bg,                       # background
                σ, σ,                             # σ_x, σ_y
                emitter.σ_photons, emitter.σ_bg;  # σ_photons, σ_bg
                frame=emitter.frame,
                dataset=emitter.dataset,
                track_id=emitter.track_id,
                id=emitter.id
            )
        end
    end
    
    metadata = copy(smld.metadata)
    metadata["simulation_type"] = "noisy_model"
    metadata["psf_width"] = σ_psf
    
    return BasicSMLD(
        new_emitters,
        smld.camera,
        smld.n_frames,
        smld.n_datasets,
        metadata
    )
end
