# test_README.jl - Executable mirror of the code examples in README.md.
#
# Run this to verify the README's examples work against the current API
# (v0.7.0: *Config parameter objects + 2-tuple `simulate`/`gen_images` returns).
# Uses CairoMakie (headless-friendly) for the visualization example.

using Pkg
Pkg.activate(@__DIR__)
using Revise

using SMLMSim
using MicroscopePSFs
using CairoMakie
using Statistics

#==========================================================================
Quick Start: Static SMLM Simulation
==========================================================================#

camera = IdealCamera(128, 128, 0.1)                   # 128×128 pixels, 100nm pixels
params = StaticSMLMConfig(density=1.0, σ_psf=0.13)    # density 1/μm², PSF 130nm

smld_noisy, info = simulate(
    params;                                           # ';' separates positional from keyword args
    pattern=Nmer2D(n=8, d=0.1),                       # 100nm diameter ring
    camera=camera
)
# info.smld_true and info.smld_model hold the intermediate results
println("Static: generated $(length(smld_noisy.emitters)) localizations.")

#==========================================================================
Quick Start: Diffusion & Interaction Simulation
==========================================================================#

diff_params = DiffusionSMLMConfig(
    density = 0.5,           # molecules per μm²
    box_size = 10.0,         # μm
    diff_monomer = 0.1,      # μm²/s
    k_off = 0.2,             # s⁻¹ dimer dissociation rate
    dt = 0.01,               # s simulation timestep
    t_max = 10.0,            # s total simulation time
    camera_framerate = 10.0, # 10 fps (100ms per frame)
    camera_exposure = 0.1    # 100ms exposure integrates 10 timesteps per frame
)
smld_diff, diff_info = simulate(diff_params)          # returns (BasicSMLD, SimInfo)
println("Diffusion: simulated $(diff_params.t_max) s -> $(smld_diff.n_frames) frames.")

#==========================================================================
Patterns
==========================================================================#

nmer = Nmer2D(n=8, d=0.1)                             # 8 molecules in a 100nm circle
line = Line3D(λ=5.0, endpoints=[(-1.0, 0.0, -0.5), (1.0, 0.0, 0.5)])  # 5 mols/μm

#==========================================================================
Labeling
==========================================================================#

labeling = FixedLabeling()                            # exactly 1 fluorophore per site
labeling = PoissonLabeling(1.5)                       # Poisson, avg 1.5 per site
labeling = BinomialLabeling(4, 0.8)                   # 4 attachment points, 80% each
labeling = PoissonLabeling(1.5; efficiency=0.9)       # 90% of sites get labeled

smld_lab, _ = simulate(params; pattern=Nmer2D(), labeling=PoissonLabeling(1.5))
println("Labeling: $(length(smld_lab.emitters)) localizations with Poisson labeling.")

#==========================================================================
Molecules & Photophysics
==========================================================================#

# Two-state blinking model using the positional rate-matrix constructor
fluor = GenericFluor(10000.0, [-10.0 10.0; 1e-2 -1e-2])  # γ=1e4, k_off=10, k_on=1e-2

#==========================================================================
Image Generation
==========================================================================#

psf = GaussianPSF(0.15)                               # 150nm PSF width
images, img_info = gen_images(smld_diff, psf;
    support=1.0,                                      # PSF support radius in μm (faster)
    poisson_noise=true                                # add shot noise
)
println("Image generation: $(img_info.frames_generated) images, stack size $(size(images)).")

#==========================================================================
sCMOS Camera with Realistic Noise
==========================================================================#

# sCMOS camera (128×128 pixels, 100nm pixels, 1.6 e⁻ read noise)
camera_scmos = SCMOSCamera(128, 128, 0.1, 1.6)
smld_scmos, _ = simulate(params;
    pattern=Nmer2D(n=8, d=0.1),
    camera=camera_scmos
)
# Full sCMOS noise model: QE, Poisson, read noise, gain, offset
images_scmos, _ = gen_images(smld_scmos, GaussianPSF(0.15); bg=10.0, camera_noise=true)
println("sCMOS static images: size $(size(images_scmos)).")

# For diffusion simulations
diff_scmos_params = DiffusionSMLMConfig(density=0.5, box_size=10.0)
smld_diff_scmos, _ = simulate(diff_scmos_params; camera=camera_scmos, override_count=10)
println("sCMOS diffusion: $(length(smld_diff_scmos.emitters)) emitters across all frames.")

#==========================================================================
Example Workflow: Static Simulation & Visualization
==========================================================================#

camera_viz = IdealCamera(128, 128, 0.1)
params_viz = StaticSMLMConfig(density=1.0, σ_psf=0.13)
smld_viz, _ = simulate(
    params_viz;
    pattern=Nmer2D(n=6, d=0.2),                       # hexamer
    camera=camera_viz
)

emitters = smld_viz.emitters
x_coords = [e.x for e in emitters]
y_coords = [e.y for e in emitters]
photons  = [e.photons for e in emitters]

fig = Figure(size=(600, 500))
ax = Axis(fig[1, 1],
    title="Simulated SMLM Localizations (Hexamer)",
    xlabel="x (μm)", ylabel="y (μm)",
    aspect=DataAspect(), yreversed=true
)
scatter!(ax, x_coords, y_coords, color=photons, colormap=:viridis, markersize=4, alpha=0.7)
Colorbar(fig[1, 2], colormap=:viridis, label="Photons")
display(fig)
# save("smlm_hexamer.png", fig)

#==========================================================================
Example Workflow: Diffusion with Realistic sCMOS Noise
==========================================================================#

camera_scmos2 = SCMOSCamera(64, 64, 0.1, 1.6)         # 64×64 pixels, 100nm/px, 1.6 e⁻
params_spt = DiffusionSMLMConfig(
    density = 1.0,            # 1 molecule/μm²
    box_size = 6.4,           # 6.4×6.4 μm field
    diff_monomer = 0.1,       # 0.1 μm²/s diffusion
    t_max = 0.5,              # 0.5 second total
    camera_framerate = 100.0  # 100 fps
)
smld_spt, _ = simulate(params_spt; camera=camera_scmos2, photons=200.0)

# sCMOS images vs ideal-camera images for the same data
images_spt_scmos, _ = gen_images(smld_spt, GaussianPSF(0.13); bg=10.0, camera_noise=true)

camera_ideal = IdealCamera(64, 64, 0.1)
smld_spt_ideal = BasicSMLD(smld_spt.emitters, camera_ideal, smld_spt.n_frames, smld_spt.n_datasets)
images_spt_ideal, _ = gen_images(smld_spt_ideal, GaussianPSF(0.13); bg=10.0, poisson_noise=true)

println("sCMOS: mean=$(round(mean(images_spt_scmos), digits=1)) ADU, std=$(round(std(images_spt_scmos), digits=1))")
println("Ideal: mean=$(round(mean(images_spt_ideal), digits=1)) photons, std=$(round(std(images_spt_ideal), digits=1))")

println("\nAll README examples executed successfully.")
