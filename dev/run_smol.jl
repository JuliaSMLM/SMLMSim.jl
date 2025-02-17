using Revise
using SMLMSim
using MicroscopePSFs
using CairoMakie

# 1. Basic Diffusion Simulation
println("Running basic diffusion simulation...")

# Set up simulation parameters
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

# Basic visualization
visualize_sequence(systems, filename="diffusion.mp4", framerate=30)

# 2. Generate Microscope Images
println("Generating microscope images...")

# Create PSF model
psf = Gaussian2D(0.15)  # σ = 150nm

# Generate full sequence of images
images = SMLMSim.InteractionDiffusion.gen_image_sequence(psf, systems, frame_integration=10)

# Generate dimer-only images
dimer_systems = get_dimers(systems)
dimer_images = SMLMSim.InteractionDiffusion.gen_image_sequence(psf, dimer_systems, frame_integration=10)

# Visualize example frames
fig = Figure(resolution=(800, 400))

# Show full system
ax1 = Axis(fig[1, 1], title="All molecules")
heatmap!(ax1, images[:, :, 1], colormap=:viridis)

# Show dimers only
ax2 = Axis(fig[1, 2], title="Dimers only")
heatmap!(ax2, dimer_images[:, :, 1], colormap=:viridis)

save("comparison_frame.png", fig)

# 3. Analysis
println("Analyzing dimer dynamics...")

# Calculate dimer fraction over time
dimer_fractions = [count(m -> m.state == 2, sys.molecules) / length(sys.molecules) 
                   for sys in systems]
times = range(0, params.t_max, length=length(systems))

# Plot dimer fraction
fig2 = Figure()
ax = Axis(fig2[1,1], 
    title="Dimer Fraction Over Time",
    xlabel="Time (s)",
    ylabel="Fraction of Molecules in Dimers")
lines!(ax, times, dimer_fractions)
save("dimer_fraction.png", fig2)

println("Demo complete! Check output files:
- diffusion.mp4: Full simulation visualization
- comparison_frame.png: Example microscope images
- dimer_fraction.png: Dimer fraction analysis")

# Example of accessing individual molecule data
println("\nExample molecular data from final frame:")
final_system = systems[end]
for (i, mol) in enumerate(final_system.molecules[1:3])  # First 3 molecules
    println("Molecule $i:")
    println("  Position: ($(round(mol.x, digits=3)), $(round(mol.y, digits=3)))")
    println("  State: $(mol.state == 1 ? "Monomer" : "Dimer")")
    if mol.state == 2
        println("  Linked to molecule: $(mol.link)")
    end
end