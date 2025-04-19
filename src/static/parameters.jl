"""
    StaticSMLMParams <: SMLMSimParams

Parameters for static SMLM simulation.

# Fields
- `density::Float64`: density in particles per square micron
- `σ_psf::Float64`: PSF width in microns
- `minphotons::Int`: minimum photons for detection
- `ndatasets::Int`: number of datasets to simulate
- `nframes::Int`: number of frames per dataset
- `framerate::Float64`: frames per second
- `ndims::Int`: dimensionality (2 or 3)
- `zrange::Vector{Float64}`: axial range for 3D simulations [min_z, max_z]

# Examples
```julia
# Default parameters
params = StaticSMLMParams()

# Custom parameters
params = StaticSMLMParams(
    density = 2.0,              # 2 particles per μm²
    σ_psf = 0.15,         # 150nm PSF width
    minphotons = 100,     # minimum photons for detection
    ndatasets = 5,        # 5 independent datasets
    nframes = 2000,       # 2000 frames per dataset
    framerate = 100.0,    # 100 frames per second
    ndims = 3,            # 3D simulation
    zrange = [-2.0, 2.0]  # 4μm axial range
)
```
"""
Base.@kwdef mutable struct StaticSMLMParams <: SMLMSimParams
    density::Float64 = 1.0
    σ_psf::Float64 = 0.13
    minphotons::Int = 50
    ndatasets::Int = 10
    nframes::Int = 1000
    framerate::Float64 = 50.0
    ndims::Int = 2
    zrange::Vector{Float64} = [-1.0, 1.0]
    
    function StaticSMLMParams(
        density, σ_psf, minphotons, ndatasets, nframes, framerate, ndims, zrange
    )
        # Input validation
        if density <= 0
            throw(ArgumentError("Density must be positive"))
        end
        
        if σ_psf <= 0
            throw(ArgumentError("PSF width (σ_psf) must be positive"))
        end
        
        if minphotons < 0
            throw(ArgumentError("Minimum photon count must be non-negative"))
        end
        
        if ndatasets <= 0
            throw(ArgumentError("Number of datasets must be positive"))
        end
        
        if nframes <= 0
            throw(ArgumentError("Number of frames must be positive"))
        end
        
        if framerate <= 0
            throw(ArgumentError("Frame rate must be positive"))
        end
        
        if ndims != 2 && ndims != 3
            throw(ArgumentError("Number of dimensions must be 2 or 3"))
        end
        
        if length(zrange) != 2 || zrange[1] >= zrange[2]
            throw(ArgumentError("zrange must be a vector of two values [min_z, max_z] where min_z < max_z"))
        end
        
        new(density, σ_psf, minphotons, ndatasets, nframes, framerate, ndims, zrange)
    end
end

function Base.show(io::IO, params::StaticSMLMParams)
    println(io, "StaticSMLMParams:")
    println(io, "  density    = $(params.density) particles/μm²")
    println(io, "  σ_psf     = $(params.σ_psf) μm")
    println(io, "  minphotons = $(params.minphotons)")
    println(io, "  ndatasets = $(params.ndatasets)")
    println(io, "  nframes   = $(params.nframes)")
    println(io, "  framerate = $(params.framerate) Hz")
    println(io, "  ndims     = $(params.ndims)D")
    print(io,   "  zrange    = [$(params.zrange[1]), $(params.zrange[2])] μm")
end

function Base.show(io::IO, ::MIME"text/plain", params::StaticSMLMParams)
    println(io, "StaticSMLMParams:")
    println(io, "  Density           = $(params.density) particles/μm²")
    println(io, "  PSF width (σ)     = $(params.σ_psf) μm ($(params.σ_psf*1000) nm)")
    println(io, "  Min photons       = $(params.minphotons)")
    println(io, "  Datasets          = $(params.ndatasets)")
    println(io, "  Frames per dataset= $(params.nframes)")
    println(io, "  Frame rate        = $(params.framerate) Hz ($(1000/params.framerate) ms/frame)")
    println(io, "  Dimensions        = $(params.ndims)D")
    print(io,   "  Z-range           = [$(params.zrange[1]), $(params.zrange[2])] μm")
end