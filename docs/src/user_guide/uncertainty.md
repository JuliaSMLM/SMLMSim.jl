```@meta
CurrentModule = SMLMSim
```

# Localization Uncertainty

This page explains how localization uncertainty is implemented in SMLMSim and how to customize it for realistic SMLM simulations.

## Overview

In SMLM, localization uncertainty arises from:

1. **Photon statistics**: Finite photon counts lead to statistical errors in position estimation
2. **Background noise**: Reduces precision in determining the true position
3. **Pixelation**: Camera pixel size impacts the precision of measurements
4. **PSF model mismatch**: Differences between actual and fitted PSF models

SMLMSim primarily models the photon statistics component, which is typically the dominant source of uncertainty in SMLM.

## Theoretical Basis

The fundamental limit of localization precision is given by the Cramér-Rao Lower Bound (CRLB). For a Gaussian PSF in the presence of background noise, the uncertainty in each dimension can be approximated as:

$$\sigma_{\text{loc}} \approx \frac{\sigma_{\text{PSF}}}{\sqrt{N}}$$

Where:
- $\sigma_{\text{PSF}}$ is the width of the point spread function
- $N$ is the number of photons detected from the emitter

This formula captures the key insight that localization precision improves with:
- More photons (higher signal)
- Smaller PSF width

## Implementing Uncertainty in SMLMSim

SMLMSim applies localization uncertainty using the `noise()` function, which adds position errors based on photon counts:

```julia
# Apply localization uncertainty to model with kinetic blinking
smld_noisy = noise(smld_model, 0.13)  # 0.13μm = 130nm PSF width
```

For 3D simulations, you specify the PSF width in each dimension:

```julia
# Apply 3D localization uncertainty
smld_noisy_3d = noise(smld_model_3d, [0.13, 0.13, 0.39])  # [σx, σy, σz] in μm
```

Note that the axial (z) uncertainty is typically 2-3× larger than the lateral (x,y) uncertainty.

## Uncertainty in Simulations

The high-level `simulate()` function handles uncertainty automatically:

```julia
smld_true, smld_model, smld_noisy = simulate(
    σ_psf=0.13,           # PSF width in μm (130nm)
    pattern=Nmer2D(),
    camera=IdealCamera(128, 128, 0.1)
)
```

The third returned object (`smld_noisy`) contains emitters with:
- Position noise added according to the PSF width and photon counts
- Uncertainty fields (`σ_x`, `σ_y`, and for 3D `σ_z`) populated with the theoretical uncertainty values

## Accessing Uncertainty Values

You can access the uncertainty values for each emitter:

```julia
# Get position uncertainties for all emitters
σ_x_values = [e.σ_x for e in smld_noisy.emitters]
σ_y_values = [e.σ_y for e in smld_noisy.emitters]

# For 3D data
σ_z_values = [e.σ_z for e in smld_noisy.emitters]

# Plot relationship between photons and uncertainty
using CairoMakie
photons = [e.photons for e in smld_noisy.emitters]
fig = Figure(size=(600, 400))
ax = Axis(fig[1, 1], 
    xlabel="Photons", ylabel="σx (μm)",
    xscale=log10, yscale=log10,
    title="Localization Uncertainty")
scatter!(ax, photons, σ_x_values)
fig
```

## Customizing Uncertainty Models

### PSF Width

The PSF width is the most important parameter affecting localization uncertainty:

```julia
# Simulation with wider PSF (150nm)
smld_true, smld_model, smld_wide_psf = simulate(σ_psf=0.15)

# Simulation with narrower PSF (100nm)
smld_true, smld_model, smld_narrow_psf = simulate(σ_psf=0.10)
```

Realistic PSF widths depend on:
- Wavelength (λ)
- Numerical aperture (NA)
- Optical aberrations

For visible light microscopy, typical values range from 100-250nm.

### Photon Counts

The number of photons is controlled by:

1. **Fluorophore emission rate**: Set in the molecule model
2. **Exposure time**: Determined by the framerate
3. **Detection threshold**: Set with minphotons

```julia
using SMLMSim

# Define bright and dim fluorophores using the positional constructor
bright_fluor = GenericFluor(5e4, [-5.0 5.0; 1.0 -1.0]) # γ=5e4, k_off=5, k_on=1
dim_fluor = GenericFluor(5e3, [-5.0 5.0; 1.0 -1.0])   # γ=5e3, k_off=5, k_on=1

# Bright emitters with lower uncertainty
smld_true, smld_model, smld_bright = simulate(molecule=bright_fluor)

# Dim emitters with higher uncertainty
smld_true, smld_model, smld_dim = simulate(molecule=dim_fluor)
```

### Advanced: Custom Uncertainty Models

For more complex uncertainty models, you can modify the `noise()` function directly:

```julia
# Example: Add non-Gaussian errors
function custom_noise(smld, σ_psf)
    # Start with standard Gaussian noise
    smld_noisy = noise(smld, σ_psf)
    
    # Add occasional large errors (e.g., fitting failures)
    for (i, e) in enumerate(smld_noisy.emitters)
        if rand() < 0.05  # 5% chance of large error
            # Add additional displacement
            smld_noisy.emitters[i].x += randn() * 0.1  # 100nm additional error
            smld_noisy.emitters[i].y += randn() * 0.1
        end
    end
    
    return smld_noisy
end

# Use custom noise model
smld_custom = custom_noise(smld_model, 0.13)
```

## Visualizing Localization Uncertainty

To visualize the effect of localization uncertainty:

```julia
using CairoMakie

# Compare ground truth with noisy localizations
function plot_comparison(smld_true, smld_noisy)
    # Extract coordinates
    x_true = [e.x for e in smld_true.emitters]
    y_true = [e.y for e in smld_true.emitters]
    
    x_noisy = [e.x for e in smld_noisy.emitters]
    y_noisy = [e.y for e in smld_noisy.emitters]
    
    # Create figure and plot
    fig = Figure(size=(700, 500))
    ax = Axis(fig[1, 1], 
        xlabel="x (μm)", 
        ylabel="y (μm)", 
        aspect=DataAspect(),
        title="Localization Uncertainty")
    
    # Plot ground truth and noisy localizations
    scatter!(ax, x_true, y_true, 
        label="Ground Truth", 
        markersize=8, 
        color=:black, 
        alpha=0.7)
    
    scatter!(ax, x_noisy, y_noisy, 
        label="Noisy Localizations", 
        markersize=4, 
        color=:red, 
        alpha=0.3)
    
    axislegend(ax)
    
    return fig
end

# Use with simulation results
plot_comparison(smld_true, smld_noisy)
```

## Realistic Uncertainty Considerations

When designing realistic simulations, consider:

1. **Depth-dependent uncertainty**: In 3D imaging, uncertainty typically increases with distance from the focal plane
2. **Fluorophore orientation effects**: Dipole orientation can affect PSF shape and localization precision
3. **Sample-dependent background**: Higher background in dense regions increases uncertainty
4. **Drift**: Sample drift adds systematic errors over time

Advanced simulations could incorporate these effects for more accurate modeling of experimental conditions.