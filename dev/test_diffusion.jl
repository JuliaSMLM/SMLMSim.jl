# test_diffusion.jl - Simulate diffusing and interacting molecules

using Pkg
Pkg.activate("dev")
using Revise

using SMLMSim
using SMLMSim.InteractionDiffusion
using MicroscopePSFs
using GLMakie
using Statistics
include("stack_viewer.jl")

println("Setting up diffusion simulation parameters...")

# Create simulation parameters
temporal_sampling = 10
frame_rate = 100.0
params = DiffusionSMLMParams(
    density = 1.0,        # molecules per μm²
    box_size = 5.0,       # 5μm box 
    diff_monomer = 0.1,   # Diffusion coefficient for monomers (μm²/s)
    diff_dimer = 0.05,    # Diffusion coefficient for dimers (μm²/s)
    k_off = 1.0,          # Dissociation rate (s⁻¹)
    r_react = 0.015,      # Reaction radius (μm)
    d_dimer = 0.03,       # Dimer separation distance (μm)
    dt = 1/frame_rate/temporal_sampling,  # Time step (s)
    t_max = 5.0,          # Total simulation time (s) 
    boundary = "periodic", # Boundary condition
    ndims = 2,             # 2D simulation
    camera_framerate = frame_rate, # Camera framerate (Hz) - for generating images
)

println("Simulation parameters:")
println("- Density: $(params.density) molecules/μm²")
println("- Box size: $(params.box_size)×$(params.box_size) μm²")
println("- Expected molecules: $(round(Int, params.density * params.box_size^2))")
println("- Diffusion coefficients: D_monomer=$(params.diff_monomer) μm²/s, D_dimer=$(params.diff_dimer) μm²/s")
println("- Reaction parameters: r_react=$(params.r_react) μm, k_off=$(params.k_off) s⁻¹")
println("- Time parameters: dt=$(params.dt) s, t_max=$(params.t_max) s")

println("Running diffusion simulation...")
systems = simulate(params)
n_frames = length(systems)

println("Simulation complete with $(n_frames) frames")

# Create a PSF model for generating camera images
println("Creating PSF model...")
psf = GaussianPSF(0.13)  # 130nm PSF width

# Generate camera images
println("Generating camera images...")
images = gen_images(systems, psf;
    support=1.0,  # PSF support radius in μm
    bg=2.0,             # Background photons
    poisson_noise=true   # Add realistic Poisson noise
)

println("Image stack generated with size: $(size(images))")

# Launch the stack viewer
println("Launching stack viewer...")
fig1 = view_stack(images, title="SMLM 2D Simulation (Hexamer Pattern)")

# Show dimers only 
system_dimers = get_dimers(systems)
images_dimers = gen_images(system_dimers, psf;
    support=1.0,  # PSF support radius in μm
    bg=2.0,             # Background photons
    poisson_noise=true   # Add realistic Poisson noise
)
fig2 = view_stack(images_dimers, title="SMLM 2D Simulation (Hexamer Pattern)")

# Show monomers only
system_monomers = get_monomers(systems)
images_monomers = gen_images(system_monomers, psf;
    support=1.0,  # PSF support radius in μm
    bg=2.0,             # Background photons
    poisson_noise=true   # Add realistic Poisson noise
)
fig3 = view_stack(images_monomers, title="SMLM 2D Simulation (Hexamer Pattern)")

dimer_lifetime = analyze_dimer_lifetime(systems)
println("Dimer lifetime: $(dimer_lifetime)")


println("Analysis complete.")
