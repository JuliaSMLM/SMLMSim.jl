"""
This module provides simulation tools for diffusion
    and interaction between particles.
"""
module InteractionDiffusion

using Distributions
using CairoMakie
using MicroscopePSFs
using SMLMSim

include("types.jl")
include("diffusion.jl")
include("smoluchowski.jl")
include("visualize.jl")
include("microscope.jl")
include("dimer.jl")

export smoluchowski, gen_movie, show_frame, gen_image, gen_image_stack

end
