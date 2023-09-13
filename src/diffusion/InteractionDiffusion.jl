"""
This module provides simulation tools for diffusion
    and interaction between particles.
"""
module InteractionDiffusion

using Distributions

include("diffusion.jl")
include("smoluchowski.jl")

export smoluchowski

end
