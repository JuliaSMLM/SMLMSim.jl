"""
    Stepper

Closed-loop simulation engine: a [`SimWorld`](@ref) holds the emitters of one or more
[`Population`](@ref)s and [`step!`](@ref) advances it over one camera exposure, returning
the noise-free expected photons per pixel. Camera noise is a separate call
([`scmos_noise!`](@ref) or [`poisson_noise!`](@ref)) on the caller's own RNG.

# Usage
```julia
using SMLMSim.Stepper
```
"""
module Stepper

using Random
using LinearAlgebra
using SMLMData: AbstractCamera, IdealCamera, SCMOSCamera
using MicroscopePSFs: GaussianPSF
using ..Core: GenericFluor
using ..CameraImages: RenderBuffer, render_gaussian!, StampTable, render_stamp!

include("types.jl")
include("kinetics.jl")
include("background.jl")
include("step.jl")

export Population, BackgroundModel, SimWorld, UniformExcitation
export step!, layers, next_switch, gen_background

end # module
