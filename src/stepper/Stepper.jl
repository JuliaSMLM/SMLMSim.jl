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
using SMLMData: AbstractCamera, IdealCamera, SCMOSCamera
using MicroscopePSFs: GaussianPSF
using ..Core: GenericFluor
using ..InteractionDiffusion: DiffusionSMLMConfig
using ..CameraImages: _uniform_pitch, RenderBuffer, render_gaussian!, StampTable, render_stamp!

include("types.jl")
include("kinetics.jl")
include("dimers.jl")
include("background.jl")
include("truth.jl")
include("step.jl")

export Population, DimerKinetics, BackgroundModel, SimWorld, UniformExcitation, EvanescentExcitation
export FrameTruth, frame_truth, next_switch, gen_background  # step! and layers are public, not exported (generic names)

end # module
