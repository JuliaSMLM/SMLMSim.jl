# using Pkg
# Pkg.activate("dev")

# using Revise
using SMLMSim # must be 'dev'ed in the current environment
using CairoMakie

#==========================================================================
Basic Usage
==========================================================================#

# Basic simulation with default parameters
camera = IdealCamera(1:128, 1:128, 0.1)  # 128×128 pixels, 100nm pixels
smld_true, smld_model, smld_noisy = simulate(
    camera=camera
)

# More customized simulation
smld_true, smld_model, smld_noisy = simulate(;
    ρ=1.0,                # emitters per μm²
    σ_psf=0.13,           # PSF width in μm (130nm)
    minphotons=50,        # minimum photons for detection
    ndatasets=10,         # number of independent datasets
    nframes=1000,         # frames per dataset
    framerate=50.0,       # frames per second
    pattern=Nmer2D(n=6, d=0.2),  # hexamer with 200nm diameter
    molecule=GenericFluor(; q=[0 50; 1e-2 0]),  # rates in 1/s
    camera=IdealCamera(; ypixels=256, xpixels=128, pixelsize=0.1)  # pixelsize in μm
)

#==========================================================================
Pattern Types
==========================================================================#

# 2D Patterns
# N molecules arranged in a circle
nmer = Nmer2D(n=8, d=0.1)  # 8 molecules in a 100nm diameter circle

# Linear pattern with random positions
line = Line2D(λ=5.0, endpoints=[(-2.0, 0.0), (2.0, 0.0)])  # 5 molecules per μm along line

# 3D Patterns
# N molecules arranged in a circle at z=0
nmer3d = Nmer3D(n=8, d=0.1)  # 8 molecules in a 100nm diameter circle

# 3D line with random positions
line3d = Line3D(λ=5.0, endpoints=[(-1.0, 0.0, -0.5), (1.0, 0.0, 0.5)])

#==========================================================================
Molecule Models
==========================================================================#

# Generic fluorophore with two-state kinetics
fluor = GenericFluor(
    γ=10000.0,           # photon emission rate in Hz
    q=[0 10; 1e-2 0]     # transition rate matrix: state 1 ↔ state 2
)

#==========================================================================
Diffusion and Interaction Simulation
==========================================================================#

# Set up parameters for Smoluchowski diffusion simulation
params = SmoluchowskiParams(
    density = 0.5,        # molecules per μm²
    box_size = 10.0,      # μm
    diff_monomer = 0.1,   # μm²/s
    diff_dimer = 0.05,    # μm²/s
    k_off = 0.2,          # s⁻¹
    r_react = 0.01,       # μm
    d_dimer = 0.05,       # μm
    dt = 0.01,            # s
    t_max = 10.0          # s
)

# Run simulation
systems = simulate(params)

# Visualize the simulation
visualize_sequence(systems, filename="diffusion.mp4", framerate=30)

# Generate microscope images
psf = Gaussian2D(0.15)  # 150nm PSF width
images = gen_image_sequence(
    psf, 
    systems,
    frame_integration=10
)

# Extract only dimers
dimer_systems = get_dimers(systems)
dimer_images = gen_image_sequence(
    psf, 
    dimer_systems, 
    frame_integration=10
)

#==========================================================================
2D Simulation with Visualization
==========================================================================#

# Create camera with physical pixel size
camera_viz = IdealCamera(1:128, 1:256, 0.1)  # 128×256 pixels, 100nm pixels

# Simulation parameters in physical units
smld_true_viz, smld_model_viz, smld_noisy_viz = simulate(;
    ρ=1.0,                # emitters per μm²
    σ_psf=0.13,           # PSF width in μm
    pattern=Nmer2D(n=6, d=0.2),  # hexamer with 200nm diameter
    camera=camera_viz
)

# Extract coordinates from emitters
x_noisy = [e.x for e in smld_noisy_viz.emitters]
y_noisy = [e.y for e in smld_noisy_viz.emitters]
photons = [e.photons for e in smld_noisy_viz.emitters]

# Create figure and plot results
fig = Figure(size=(800, 600))
ax = Axis(fig[1, 1], 
    title="Simulated SMLM Localizations",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect()
)

# Scatter plot with photon counts as color
scatter!(ax, x_noisy, y_noisy, 
    color=photons,
    colormap=:viridis,
    markersize=4,
    alpha=0.6
)

Colorbar(fig[1, 2], colormap=:viridis, label="Photons")

# Show or save the figure
display(fig)
# save("smlm_simulation.png", fig)

#==========================================================================
3D Simulation
==========================================================================#

# Create camera with physical pixel size
camera_3d = IdealCamera(1:128, 1:256, 0.1)  # 128×256 pixels, 100nm pixels

# Simulation parameters in physical units
smld_true_3d, smld_model_3d, smld_noisy_3d = simulate(;
    ρ=0.5,                # emitters per μm²
    pattern=Nmer3D(n=8, d=0.3),  # 3D pattern
    camera=camera_3d,
    zrange=[-2.0, 2.0]    # 4μm axial range
)