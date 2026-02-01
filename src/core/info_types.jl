"""
Info types for simulation results following the (result, info) return convention.

These structs capture metadata about simulation execution including timing,
counts, and intermediate results.
"""

"""
    SimInfo

Metadata about a simulation run, returned as the second element of the
(result, info) tuple from `simulate()`.

# Fields
## Common fields (ecosystem convention)
- `elapsed_ns::UInt64`: Execution time in nanoseconds
- `backend::Symbol`: Computation backend (`:cpu`)
- `device_id::Int`: Device identifier (-1 for CPU)

## Simulation-specific
- `seed::Union{UInt64, Nothing}`: RNG seed for reproducibility (if set)

## Static simulation outputs
- `smld_true::Union{SMLD, Nothing}`: Ground truth positions
- `smld_model::Union{SMLD, Nothing}`: After kinetic model (blinking)

## Counts
- `n_patterns::Int`: Number of spatial patterns generated
- `n_emitters::Int`: Total emitters in ground truth
- `n_localizations::Int`: Total localizations after kinetic model
- `n_frames::Int`: Number of frames simulated
"""
struct SimInfo
    # Common fields (ecosystem convention)
    elapsed_ns::UInt64
    backend::Symbol
    device_id::Int

    # Simulation-specific
    seed::Union{UInt64, Nothing}

    # Static simulation outputs (using abstract SMLD type)
    smld_true::Union{SMLD, Nothing}
    smld_model::Union{SMLD, Nothing}

    # Counts
    n_patterns::Int
    n_emitters::Int
    n_localizations::Int
    n_frames::Int
end

"""
    SimInfo(; elapsed_ns, smld_true=nothing, smld_model=nothing, seed=nothing,
            n_patterns=0, n_emitters=0, n_localizations=0, n_frames=0)

Convenience constructor for SimInfo with keyword arguments.
"""
function SimInfo(;
    elapsed_ns::UInt64,
    smld_true::Union{SMLD, Nothing}=nothing,
    smld_model::Union{SMLD, Nothing}=nothing,
    seed::Union{UInt64, Nothing}=nothing,
    n_patterns::Int=0,
    n_emitters::Int=0,
    n_localizations::Int=0,
    n_frames::Int=0
)
    SimInfo(
        elapsed_ns,
        :cpu,
        -1,
        seed,
        smld_true,
        smld_model,
        n_patterns,
        n_emitters,
        n_localizations,
        n_frames
    )
end

"""
    ImageInfo

Metadata about image generation, returned as the second element of the
(result, info) tuple from `gen_images()`.

# Fields
## Common fields (ecosystem convention)
- `elapsed_ns::UInt64`: Execution time in nanoseconds
- `backend::Symbol`: Computation backend (`:cpu`)
- `device_id::Int`: Device identifier (-1 for CPU)

## Image generation specific
- `frames_generated::Int`: Number of frames generated
- `n_photons_total::Float64`: Total photons across all frames
- `output_size::Tuple{Int,Int,Int}`: Output dimensions (H, W, T)
"""
struct ImageInfo
    # Common fields
    elapsed_ns::UInt64
    backend::Symbol
    device_id::Int

    # Image generation specific
    frames_generated::Int
    n_photons_total::Float64
    output_size::Tuple{Int,Int,Int}
end

"""
    ImageInfo(; elapsed_ns, frames_generated, n_photons_total, output_size)

Convenience constructor for ImageInfo with keyword arguments.
"""
function ImageInfo(;
    elapsed_ns::UInt64,
    frames_generated::Int,
    n_photons_total::Float64,
    output_size::Tuple{Int,Int,Int}
)
    ImageInfo(
        elapsed_ns,
        :cpu,
        -1,
        frames_generated,
        n_photons_total,
        output_size
    )
end

# Pretty printing
function Base.show(io::IO, info::SimInfo)
    elapsed_ms = info.elapsed_ns / 1_000_000
    print(io, "SimInfo(")
    print(io, "elapsed=$(round(elapsed_ms, digits=2))ms, ")
    print(io, "n_emitters=$(info.n_emitters), ")
    print(io, "n_localizations=$(info.n_localizations), ")
    print(io, "n_frames=$(info.n_frames)")
    print(io, ")")
end

function Base.show(io::IO, info::ImageInfo)
    elapsed_ms = info.elapsed_ns / 1_000_000
    print(io, "ImageInfo(")
    print(io, "elapsed=$(round(elapsed_ms, digits=2))ms, ")
    print(io, "frames=$(info.frames_generated), ")
    print(io, "photons=$(round(info.n_photons_total, digits=1)), ")
    print(io, "size=$(info.output_size)")
    print(io, ")")
end
