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
params = StaticSMLMConfig(density=1.0, σ_psf=0.13, nframes=1000)
camera = IdealCamera(128, 128, 0.1)  # 128×128 pixels, 100nm pixels
pattern = Nmer2D(n=6, d=0.2)  # Hexamer with 200nm diameter
molecule = GenericFluor(1e4, [-50.0 50.0; 1e-2 -1e-2])  # Blinking model

# Run simulation
smld_noisy, info = simulate(params, pattern=pattern, molecule=molecule, camera=camera)

# Create a PSF model (Gaussian with 150nm width)
psf = GaussianPSF(0.15)  # 150nm PSF width

# Generate image stack from emitter data
images, img_info = gen_images(info.smld_model, psf;
    bg=5.0,             # background photons per pixel
    poisson_noise=true  # add realistic photon counting noise
)
```

The resulting `images` is a 3D array with dimensions `[height, width, frames]` that can be used for visualization or algorithm development.

### Creating a Single Frame

To generate an image for a specific frame:

```julia
# Generate image for frame 10
frame_image, frame_info = gen_image(info.smld_model, psf, 10;
    bg=5.0,
    poisson_noise=true
)
```

## Function Parameters

The `gen_images` function has the following signature:

```julia
gen_images(smld::SMLD, psf::AbstractPSF; kwargs...) -> (Array{T, 3}, info) where T<:Real
```

### Required Parameters

- `smld::SMLD`: Single molecule localization data container
- `psf::AbstractPSF`: Point spread function model from MicroscopePSFs.jl

### Optional Keyword Arguments

```julia
# Complete example with all available options
images, img_info = gen_images(info.smld_model, psf;
    dataset::Int=1,                # Dataset number to use from SMLD
    frames=nothing,                # Specific frames to generate (default: all frames)
    support=Inf,                   # PSF support region size in μm (Inf, scalar, or tuple)
    sampling=2,                    # Supersampling factor for PSF integration
    threaded=true,                 # Enable multithreading for faster computation
    bg=0.0,                        # Background signal level (photons per pixel)
    poisson_noise=false,           # Apply Poisson noise only (simple shot noise)
    camera_noise=false             # Apply full camera noise model (requires SCMOSCamera)
)
```

## Noise Options

Realistic microscope images include various noise sources. SMLMSim supports the following options:

### Photon Shot Noise

Photon shot noise follows a Poisson distribution and is the most fundamental noise source in fluorescence microscopy:

```julia
# Enable Poisson noise (default: false)
images, img_info = gen_images(info.smld_model, psf, poisson_noise=true)
```

When enabled, each pixel's intensity is drawn from a Poisson distribution with λ equal to the expected photon count.

### Background Signal

Background noise can be added as a constant offset to all pixels:

```julia
# Add background signal (photons per pixel)
images, img_info = gen_images(info.smld_model, psf, bg=10.0)
```

When combined with Poisson noise, the background is included in the Poisson sampling process for a realistic noise model.

### sCMOS Camera Noise

For realistic camera noise modeling, SMLMSim supports sCMOS cameras with per-pixel calibration parameters. This applies the full detection chain from photons to ADU (analog-to-digital units):

```julia
using SMLMSim
using MicroscopePSFs

# Create an sCMOS camera with realistic parameters
# Parameters: width, height, pixel_size (μm), readnoise (e⁻ RMS)
camera_scmos = SCMOSCamera(128, 128, 0.1, 1.6)  # 1.6 e⁻ RMS read noise

# Run simulation with sCMOS camera
params = StaticSMLMConfig(density=1.0, σ_psf=0.13)
smld_noisy, info = simulate(
    params,
    pattern=Nmer2D(n=8, d=0.1),
    camera=camera_scmos
)

# Generate images with full sCMOS noise model
psf = GaussianPSF(0.15)
images, img_info = gen_images(info.smld_model, psf, bg=10.0, camera_noise=true)
```

The sCMOS noise model applies these transformations in order:

1. **Quantum Efficiency (QE)**: Converts photons to photoelectrons (typically 70-95%)
2. **Poisson Noise**: Shot noise on the photoelectron count
3. **Read Noise**: Gaussian noise per pixel (amplifier/ADC noise)
4. **Gain**: Converts electrons to ADU (e.g., 0.5 e⁻/ADU)
5. **Offset**: Adds dark level baseline (e.g., 100 ADU)

#### Advanced sCMOS Configuration

You can specify all calibration parameters explicitly:

```julia
camera_scmos = SCMOSCamera(
    128, 128, 0.1;
    offset=100.0,      # 100 ADU dark level
    gain=0.5,          # 0.5 e⁻/ADU conversion
    readnoise=1.6,     # 1.6 e⁻ RMS read noise
    qe=0.95            # 95% quantum efficiency
)
```

Each parameter can be either:
- A scalar value (uniform across all pixels)
- A matrix (spatially varying, per-pixel calibration)

#### IdealCamera vs SCMOSCamera

```julia
# IdealCamera: Use with poisson_noise for simple shot noise
camera_ideal = IdealCamera(128, 128, 0.1)
smld = BasicSMLD(emitters, camera_ideal, n_frames, n_datasets)
images_ideal, img_info_ideal = gen_images(smld, psf, bg=10.0, poisson_noise=true)

# SCMOSCamera: Use with camera_noise for realistic noise
camera_scmos = SCMOSCamera(128, 128, 0.1, 1.6)
smld = BasicSMLD(emitters, camera_scmos, n_frames, n_datasets)
images_scmos, img_info_scmos = gen_images(smld, psf, bg=10.0, camera_noise=true)
```

!!! warning "Camera Type Compatibility"
    The `camera_noise=true` parameter only works with `SCMOSCamera`. Using it with `IdealCamera` will show a warning and be ignored. For `IdealCamera`, use `poisson_noise=true` instead.

#### Spatially-Varying Calibration Maps

Real sCMOS sensors have per-pixel variation in offset, gain, and readnoise. You can model this with calibration maps:

```julia
using SMLMSim
using MicroscopePSFs

# Create per-pixel calibration maps (128×128 sensor)
n_pixels = 128

# Readnoise: mean 1.6 e⁻ with ±0.3 e⁻ variation
readnoise_map = 1.6 .+ 0.3 .* randn(n_pixels, n_pixels)
readnoise_map = max.(readnoise_map, 0.5)  # Ensure positive values

# Offset: mean 100 ADU with ±5 ADU variation plus gradient
offset_map = zeros(n_pixels, n_pixels)
for i in 1:n_pixels, j in 1:n_pixels
    gradient = 10.0 * (i / n_pixels)  # Vertical gradient
    noise = 5.0 * randn()
    offset_map[i, j] = 100.0 + gradient + noise
end

# Gain: mostly uniform with slight pattern
gain_map = 0.5 .+ 0.05 .* sin.(2π .* (1:n_pixels) ./ 20) * cos.(2π .* (1:n_pixels)' ./ 20)

# Create sCMOS camera with spatially-varying calibration
pixel_size = 0.1  # 100 nm
camera_scmos = SCMOSCamera(
    collect(0:pixel_size:(n_pixels*pixel_size)),
    collect(0:pixel_size:(n_pixels*pixel_size)),
    offset_map,
    gain_map,
    readnoise_map,
    0.95  # Uniform quantum efficiency (could also be a matrix)
)

# Generate images with realistic spatial artifacts
# (offset gradient, readnoise stripes, etc.)
images, img_info = gen_images(smld, psf, bg=10.0, camera_noise=true)
```

This creates realistic spatial artifacts visible in the images, including:
- **Offset gradients**: Fixed pattern across sensor
- **Readnoise variation**: Some pixels noisier than others
- **Gain non-uniformity**: Affects signal scaling

For demonstration of extreme artifacts, see `dev/scmos_quick_demo.jl` in the repository.

### PSF Support Region

The PSF support region controls how much of the PSF is calculated for each emitter, affecting both accuracy and performance. This parameter is specified in microns (μm):

```julia
# Define a small support region (faster but may truncate PSF wings)
images, img_info = gen_images(smld_model, psf, support=3.0)  # 3.0 μm radius around each emitter

# Define a rectangular support region
images, img_info = gen_images(smld_model, psf, support=(3.0, 3.0, 1.0, 1.0))  # (left, right, bottom, top) in μm

# Use infinite support (most accurate but slowest)
images, img_info = gen_images(smld_model, psf, support=Inf)
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

