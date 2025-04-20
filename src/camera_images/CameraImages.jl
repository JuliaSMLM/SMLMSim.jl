"""
    CameraImages

Module for generating simulated camera images from SMLM data.

This module provides functions to:
1. Generate ideal camera images by integrating emitter photons over a PSF.
2. Add realistic camera noise (e.g., Poisson noise).

# Usage
```julia
using SMLMSim.CameraImages
```
"""
module CameraImages

# All module imports should be at the top
using SMLMData
using MicroscopePSFs
using Distributions  # Required for Poisson noise functions

# Include all source files
include("gen_images.jl")
include("noise.jl")

# Export the image generation functions
export gen_images, gen_image

# Export noise functions
export poisson_noise, poisson_noise!

end # module