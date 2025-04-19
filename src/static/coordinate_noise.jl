"""
    Coordinate Noise Module

This module provides functions for modeling the localization coordinate noise
in SMLM, including detection uncertainty calculations based on photon statistics.
"""

"""
    add_coordinate_noise(emitter, σ)

Helper function to add position noise to an emitter with appropriate uncertainty.
Returns new coordinate values and uncertainty values.

For 2D emitters, σ is a scalar.
For 3D emitters, σ is a 3-element vector.
"""
function add_coordinate_noise(emitter::Emitter2DFit, σ::AbstractFloat)
    # Calculate uncertainties based on photon count
    σ_scaled = σ / sqrt(emitter.photons)
    
    # Add noise to coordinates
    x_noisy = emitter.x + randn() * σ_scaled
    y_noisy = emitter.y + randn() * σ_scaled
    
    return (x_noisy, y_noisy), (σ_scaled, σ_scaled)
end

function add_coordinate_noise(emitter::Emitter3DFit, σ::Vector{<:AbstractFloat})
    # Calculate uncertainties based on photon count
    σ_scaled = σ ./ sqrt(emitter.photons)
    
    # Add noise to coordinates
    x_noisy = emitter.x + randn() * σ_scaled[1]
    y_noisy = emitter.y + randn() * σ_scaled[2]
    z_noisy = emitter.z + randn() * σ_scaled[3]
    
    return (x_noisy, y_noisy, z_noisy), (σ_scaled[1], σ_scaled[2], σ_scaled[3])
end

"""
    apply_noise(smld::BasicSMLD, σ_psf::AbstractFloat)

Add localization uncertainty to 2D emitter positions based on photon counts.

# Arguments
- `smld::BasicSMLD`: Input SMLD containing 2D emitters
- `σ_psf::AbstractFloat`: PSF width in microns

# Returns
- `BasicSMLD`: New SMLD with noisy positions and updated uncertainties

# Example
```julia
# Then add localization noise with specific PSF width
smld_noisy = apply_noise(smld_model, 0.13)  # 130nm PSF width
```
"""
function apply_noise(smld::BasicSMLD, σ_psf::AbstractFloat)
    # Check if this is a 2D SMLD
    emitter_type = eltype(smld.emitters)
    if !(emitter_type <: Emitter2DFit)
        error("Cannot apply scalar σ_psf to non-2D emitter type: $(emitter_type)")
    end
    
    new_emitters = similar(smld.emitters)
    
    for (i, emitter) in enumerate(smld.emitters)
        if emitter.photons <= 0
            error("Emitter photon count must be positive")
        end
        
        # Add noise to coordinates and get uncertainty values
        coords, uncertainty = add_coordinate_noise(emitter, σ_psf)
        
        # Create new emitter with noisy positions
        new_emitters[i] = typeof(emitter)(
            coords...,                     # Noisy positions
            emitter.photons,               # Photons
            emitter.bg,                    # Background
            uncertainty...,                # Uncertainty estimates
            emitter.σ_photons, emitter.σ_bg;  # Original uncertainties
            frame=emitter.frame,
            dataset=emitter.dataset,
            track_id=emitter.track_id,
            id=emitter.id
        )
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

"""
    apply_noise(smld::BasicSMLD, σ_psf::Vector{<:AbstractFloat})

Add localization uncertainty to 3D emitter positions based on photon counts.

# Arguments
- `smld::BasicSMLD`: Input SMLD containing 3D emitters
- `σ_psf::Vector{<:AbstractFloat}`: PSF widths [σx, σy, σz] in microns

# Returns
- `BasicSMLD`: New SMLD with noisy positions and updated uncertainties

# Example
```julia
# Then add localization noise with specific PSF widths
σ_psf = [0.13, 0.13, 0.39]  # 130nm lateral, 390nm axial
smld_noisy = apply_noise(smld_model, σ_psf)
```
"""
function apply_noise(smld::BasicSMLD, σ_psf::Vector{<:AbstractFloat})
    # Check if this is a 3D SMLD
    emitter_type = eltype(smld.emitters)
    if !(emitter_type <: Emitter3DFit)
        error("Cannot apply vector σ_psf to non-3D emitter type: $(emitter_type)")
    end
    
    if length(σ_psf) != 3
        error("3D emitter type requires vector of 3 σ_psf values [σx, σy, σz]")
    end
    
    new_emitters = similar(smld.emitters)
    
    for (i, emitter) in enumerate(smld.emitters)
        if emitter.photons <= 0
            error("Emitter photon count must be positive")
        end
        
        # Add noise to coordinates and get uncertainty values
        coords, uncertainty = add_coordinate_noise(emitter, σ_psf)
        
        # Create new emitter with noisy positions
        new_emitters[i] = typeof(emitter)(
            coords...,                     # Noisy positions
            emitter.photons,               # Photons
            emitter.bg,                    # Background
            uncertainty...,                # Uncertainty estimates
            emitter.σ_photons, emitter.σ_bg;  # Original uncertainties
            frame=emitter.frame,
            dataset=emitter.dataset,
            track_id=emitter.track_id,
            id=emitter.id
        )
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
