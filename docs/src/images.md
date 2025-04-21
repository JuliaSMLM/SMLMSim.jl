```@meta
CurrentModule = SMLMSim
```

# Microscope Image Generation

This page explains how to generate realistic microscope images from simulated SMLM data and the available noise options.

## Overview

SMLMSim can convert SMLD objects containing emitter coordinates into realistic microscope image stacks using point spread function (PSF) models. This functionality enables:

- Creating synthetic microscope data for algorithm testing
- Visualizing simulated molecules
- Generating training data for deep learning models
- Testing localization algorithms with realistic inputs

## Basic Image Generation

The primary functions for creating images are:

- `gen_images`: Generate multiple frames (3D stack)
- `gen_image`: Generate a single frame (2D image)

### Creating a Full Image Stack

To generate a complete image stack from an SMLD object:

```julia
using SMLMSim
using MicroscopePSFs

# Create or load an SMLD object
# (This could be from either static or diffusion simulation)
params = StaticSMLMParams(density=1.0, σ_psf=0.13, nframes=1000)
camera = IdealCamera(128, 128, 0.1)  # 128×128 pixels, 100nm pixels
pattern = Nmer2D(n=6, d=0.2)  # Hexamer with 200nm diameter
molecule = GenericFluor(1e4, [-50.0 50.0; 1e-2 -1e-2])  # Blinking model

# Run simulation
smld_true, smld_model, smld_noisy = simulate(params, pattern=pattern, molecule=molecule, camera=camera)

# Create a PSF model (Gaussian with 150nm width)
psf = GaussianPSF(0.15)  # 150nm PSF width

# Generate image stack from emitter data
images = gen_images(smld_model, psf; 
    bg=5.0,             # background photons per pixel
    poisson_noise=true  # add realistic photon counting noise
)
```

The resulting `images` is a 3D array with dimensions `[height, width, frames]` that can be used for visualization or algorithm development.

### Creating a Single Frame

To generate an image for a specific frame:

```julia
# Generate image for frame 10
frame_image = gen_image(smld_model, psf, 10;
    bg=5.0,
    poisson_noise=true
)
```

## Function Parameters

The `gen_images` function has the following signature:

```julia
gen_images(smld::SMLD, psf::AbstractPSF; kwargs...) -> Array{T, 3} where T<:Real
```

### Required Parameters

- `smld::SMLD`: Single molecule localization data container
- `psf::AbstractPSF`: Point spread function model from MicroscopePSFs.jl

### Optional Keyword Arguments

```julia
# Complete example with all available options
images = gen_images(smld_model, psf;
    dataset::Int=1,                # Dataset number to use from SMLD
    frames=nothing,                # Specific frames to generate (default: all frames)
    support=Inf,                   # PSF support region size in μm (Inf, scalar, or tuple)
    sampling=2,                    # Supersampling factor for PSF integration
    threaded=true,                 # Enable multithreading for faster computation
    bg=0.0,                        # Background signal level (photons per pixel)
    poisson_noise=false,           # Apply Poisson noise
    camera_noise=false             # Apply camera read noise (not yet implemented)
)
```

## Noise Options

Realistic microscope images include various noise sources. SMLMSim supports the following options:

### Photon Shot Noise

Photon shot noise follows a Poisson distribution and is the most fundamental noise source in fluorescence microscopy:

```julia
# Enable Poisson noise (default: false)
images = gen_images(smld_model, psf, poisson_noise=true)
```

When enabled, each pixel's intensity is drawn from a Poisson distribution with λ equal to the expected photon count.

### Background Signal

Background noise can be added as a constant offset to all pixels:

```julia
# Add background signal (photons per pixel)
images = gen_images(smld_model, psf, bg=10.0)
```

When combined with Poisson noise, the background is included in the Poisson sampling process for a realistic noise model.

### PSF Support Region

The PSF support region controls how much of the PSF is calculated for each emitter, affecting both accuracy and performance. This parameter is specified in microns (μm):

```julia
# Define a small support region (faster but may truncate PSF wings)
images = gen_images(smld_model, psf, support=3.0)  # 3.0 μm radius around each emitter

# Define a rectangular support region
images = gen_images(smld_model, psf, support=(3.0, 3.0, 1.0, 1.0))  # (left, right, bottom, top) in μm

# Use infinite support (most accurate but slowest)
images = gen_images(smld_model, psf, support=Inf)
```

Reducing the support region can significantly improve performance for large datasets at the cost of some accuracy. For most PSFs, a support radius of 3-5 times the PSF width (σ) is sufficient (typically 0.5-1.0 μm for a standard SMLM PSF).

## Working with Image Stacks

The generated images can be visualized using various Julia plotting libraries:

```julia
using CairoMakie  # or GLMakie, etc.

# Display a single frame
fig = Figure()
ax = Axis(fig[1, 1], yreversed=true)
heatmap!(ax, images[:, :, 10]', colormap = :inferno)
display(fig)

# Create and save a movie from all frames
fig = Figure()
ax = Axis(fig[1, 1], yreversed=true)

# record will clear and redraw each frame for you
record(fig, "smlm_simulation.mp4", 1:size(images, 3); framerate = 10) do i
    heatmap!(ax, images[:, :, i]', colormap = :inferno)
end
```

