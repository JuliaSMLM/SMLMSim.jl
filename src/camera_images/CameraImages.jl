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
using SMLMData: AbstractCamera, IdealCamera, SCMOSCamera, SMLD, AbstractEmitter
using MicroscopePSFs
using MicroscopePSFs: SplinePSF, integrate_pixels!
using Distributions  # Required for Poisson noise functions
using Random
using SpecialFunctions: erfc

# Import ImageInfo from parent module
import ..ImageInfo

# Include all source files
include("render.jl")
include("gen_images.jl")
include("noise.jl")

# Export the image generation functions
export gen_images, gen_image

# Export noise functions
export poisson_noise, poisson_noise!, scmos_noise, scmos_noise!

# Export the erf renderer and out-of-focus stamps (shared with the Stepper)
export RenderBuffer, render_gaussian!, StampTable, render_stamp!

# Export ImageInfo (re-export from parent)
export ImageInfo

end # module