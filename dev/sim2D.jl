using Revise 
using SMLMSim
using SMLMData
using CairoMakie

# Create camera with physical pixel size
camera = IdealCamera(1:128, 1:256, 0.1)  # 128x256 pixels, 100nm pixels

# Simulation parameters in physical units
smld_true, smld_model, smld_noisy = SMLMSim.sim(;
    ρ=1.0,                # emitters per μm²
    σ_psf=0.13,          # PSF width in μm
    minphotons=50,       # minimum photons for detection
    ndatasets=10,        # number of independent datasets
    nframes=1000,        # frames per dataset
    framerate=50.0,      # frames per second
    pattern=SMLMSim.Nmer2D(),
    molecule=SMLMSim.GenericFluor(; q=[0 50; 1e-2 0]),  # rates in 1/s
    camera=camera
)

# Extract coordinates from emitters
x_noisy = [e.x for e in smld_noisy.emitters]
y_noisy = [e.y for e in smld_noisy.emitters]
photons = [e.photons for e in smld_noisy.emitters]

# Create multi-panel figure
fig = Figure(resolution=(800, 800))

# True positions
ax1 = Axis(fig[1, 1],
    title="True Positions",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect())

x_true = [e.x for e in smld_true.emitters]
y_true = [e.y for e in smld_true.emitters]
scatter!(ax1, x_true, y_true, color=:blue, markersize=4)

# Model with blinking
ax2 = Axis(fig[1, 2],
    title="Model with Blinking",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect())

x_model = [e.x for e in smld_model.emitters]
y_model = [e.y for e in smld_model.emitters]
scatter!(ax2, x_model, y_model, color=:green, markersize=4, alpha=0.5)

# Noisy localizations
ax3 = Axis(fig[2, 1:2],
    title="Noisy Localizations",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect())

# Color by photon count
scatter!(ax3, x_noisy, y_noisy, 
    color=photons,
    colormap=:viridis,
    markersize=4,
    alpha=0.5)
Colorbar(fig[2, 3], label="Photons")

# Link axes limits
linkaxes!(ax1, ax2, ax3)

# Add some stats to figure
supertitle = "Simulation Results\n" *
    "$(length(x_true)) true positions, " *
    "$(length(x_noisy)) localizations"
Label(fig[0, :], supertitle)

fig