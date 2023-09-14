"""
This module provides simulation tools for diffusion
    and interaction between particles.
"""
module InteractionDiffusion

using Distributions
using CairoMakie

include("diffusion.jl")
include("smoluchowski.jl")
include("visualize.jl")

export smoluchowski, gen_movie, show_frame

end
