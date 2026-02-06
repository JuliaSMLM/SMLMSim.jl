"""
    Info types for simulation and image generation results.

These structs provide metadata and intermediate results alongside primary outputs,
following the tuple-pattern convention: (primary_output, info).
"""

"""
    SimInfo

Metadata and intermediate results from simulation functions.

# Fields
- `elapsed_s::Float64`: Wall-clock time for simulation in seconds
- `backend::Symbol`: Computation backend (`:cpu`)
- `device_id::Int`: Device identifier (-1 for CPU)
- `seed::Union{UInt64, Nothing}`: RNG seed for reproducibility (if applicable)
- `smld_true::Union{BasicSMLD, Nothing}`: Ground truth positions (static simulation only)
- `smld_model::Union{BasicSMLD, Nothing}`: After kinetic model (static simulation only)
- `n_patterns::Int`: Number of spatial patterns generated
- `n_emitters::Int`: Total emitters simulated
- `n_localizations::Int`: Total localizations generated
- `n_frames::Int`: Number of frames in simulation

# Example
```julia
smld_noisy, info = simulate(params)
println("Simulation took \$(info.elapsed_s) seconds")
println("Generated \$(info.n_localizations) localizations from \$(info.n_emitters) emitters")

# Access intermediate results for analysis
if info.smld_true !== nothing
    # Compare true vs noisy positions
end
```
"""
struct SimInfo <: AbstractSMLMInfo
    # Common fields (ecosystem convention)
    elapsed_s::Float64
    backend::Symbol
    device_id::Int

    # Simulation-specific
    seed::Union{UInt64, Nothing}

    # Static simulation outputs (use Any to avoid type complexity)
    smld_true::Any  # Union{BasicSMLD, Nothing}
    smld_model::Any  # Union{BasicSMLD, Nothing}

    # Counts
    n_patterns::Int
    n_emitters::Int
    n_localizations::Int
    n_frames::Int
end

# Convenience constructor with defaults
function SimInfo(;
    elapsed_s::Float64=0.0,
    backend::Symbol=:cpu,
    device_id::Int=-1,
    seed::Union{UInt64, Nothing}=nothing,
    smld_true=nothing,
    smld_model=nothing,
    n_patterns::Int=0,
    n_emitters::Int=0,
    n_localizations::Int=0,
    n_frames::Int=0
)
    SimInfo(elapsed_s, backend, device_id, seed, smld_true, smld_model,
            n_patterns, n_emitters, n_localizations, n_frames)
end

"""
    ImageInfo

Metadata from image generation functions.

# Fields
- `elapsed_s::Float64`: Wall-clock time for image generation in seconds
- `backend::Symbol`: Computation backend (`:cpu`)
- `device_id::Int`: Device identifier (-1 for CPU)
- `frames_generated::Int`: Number of frames generated
- `n_photons_total::Float64`: Total photon count across all frames
- `output_size::Tuple{Int,Int,Int}`: Image dimensions (height, width, n_frames)

# Example
```julia
images, info = gen_images(smld, psf)
println("Generated \$(info.frames_generated) frames in \$(info.elapsed_s) seconds")
println("Total photons: \$(info.n_photons_total)")
```
"""
struct ImageInfo <: AbstractSMLMInfo
    # Common fields
    elapsed_s::Float64
    backend::Symbol
    device_id::Int

    # Image generation specific
    frames_generated::Int
    n_photons_total::Float64
    output_size::Tuple{Int,Int,Int}
end

# Convenience constructor with defaults
function ImageInfo(;
    elapsed_s::Float64=0.0,
    backend::Symbol=:cpu,
    device_id::Int=-1,
    frames_generated::Int=0,
    n_photons_total::Float64=0.0,
    output_size::Tuple{Int,Int,Int}=(0,0,0)
)
    ImageInfo(elapsed_s, backend, device_id, frames_generated, n_photons_total, output_size)
end
