# test_smlm2d.jl - Generate 2D SMLM data with Nmer pattern and view using stack viewer

using Pkg
Pkg.activate("dev")
using Revise

using SMLMSim
using MicroscopePSFs
using GLMakie
using Statistics
include("stack_viewer.jl")


println("Setting up simulation parameters...")

# Create a camera with 100nm pixels
pixelsize = 0.1  # 100nm pixels
camera = IdealCamera(1:256, 1:256, pixelsize)

# Create a hexamer pattern with 200nm diameter
pattern = Nmer2D(n=6, d=0.2)  # 6 molecules in a 200nm circle

# Create static SMLM parameters
params = StaticSMLMParams(
    ρ = 1.0,           # 1 pattern per square micron
    σ_psf = 0.13,      # 130nm PSF width (realistic for visible light)
    minphotons = 50,   # Default minimum photon count for detection
    ndatasets = 1,     # Generate 1 dataset
    nframes = 1000,    # Generate 1000 frames
    framerate = 10.0,  # 10 frames per second
    ndims = 2,         # 2D simulation
    zrange = [-1.0, 1.0]  # Not used for 2D, but required parameter
)

println("Running SMLM simulation with hexamer pattern...")

# Create a fluorophore model with custom blinking kinetics

# Simulation parameters from test_photophysics.jl
nframes = params.nframes
framerate = params.framerate
total_time = nframes / framerate  # Total time in seconds

# Kinetics parameters
k_off = framerate        # Off state rate (Hz) - corresponds to 1 frame duration
n_blinks_per_fluor = 10  # Desired number of events per 1000 frames
k_on = n_blinks_per_fluor / total_time  # On state rate (Hz)

println("Fluorophore parameters:")
println("- Off rate: $k_off Hz (average on time: $(1/k_off) seconds)")
println("- On rate: $k_on Hz (average off time: $(1/k_on) seconds)")
println("- Expected events per $nframes frames: $n_blinks_per_fluor")
fluor = GenericFluor(
    1e3,                        # Photon emission rate (Hz)
    [-k_off k_off; k_on -k_on]    # Rate matrix: [k_on→off, k_off→on; k_on←off, k_off←on]
                                
)

# Run simulation with custom pattern and fluorophore model
smld_true, smld_model, smld_noisy = simulate(
    params;
    pattern=pattern,
    molecule=fluor,
    camera=camera
)

println("SMLM simulation complete")
println("- True emitters: $(length(smld_true.emitters))")
println("- Model emitters after blinking: $(length(smld_model.emitters))")
println("- Noisy emitters with localization uncertainty: $(length(smld_noisy.emitters))")

# Plot the emitters 
fig1 = Figure()
ax1_1 = Axis(fig1[1, 1], title="SMLM Simulation (Hexamer Pattern)", 
    xlabel="X (μm)", ylabel="Y (μm)", yreversed=true, aspect=DataAspect())
    x = [smld_true.emitters[i].x for i in 1:length(smld_true.emitters)]
    y = [smld_true.emitters[i].y for i in 1:length(smld_true.emitters)]

scatter!(ax1_1, x,y, color=:blue, label="True Emitters")
# noisy emitters
ax1_2 = Axis(fig1[1, 2], title="SMLM Simulation Noisy Emitters", 
    xlabel="X (μm)", ylabel="Y (μm)", yreversed=true, aspect=DataAspect())
x_noisy = [smld_noisy.emitters[i].x for i in 1:length(smld_noisy.emitters)]
y_noisy = [smld_noisy.emitters[i].y for i in 1:length(smld_noisy.emitters)]
    
scatter!(ax1_2, x_noisy,y_noisy, color=:red, label="Noisy Emitters")

# link axes
linkaxes!(ax1_1, ax1_2)
display(fig1)

# Create a PSF model for generating camera images
println("Creating PSF model...")
psf = MicroscopePSFs.GaussianPSF(0.13)  # 130nm PSF width

# Generate camera images
println("Generating camera images...")
images = gen_images(smld_noisy, psf;
    support = 1.0,  # PSF support radius in μm
    bg=2.0,             # Background photons
    poisson_noise=true   # Add realistic Poisson noise
)

println("Image stack generated with size: $(size(images))")

# Launch the stack viewer
println("Launching stack viewer...")
fig2 = view_stack(images, title="SMLM 2D Simulation (Hexamer Pattern)")

println("Stack viewer is now interactive. Close the window to exit.")

# Count fluorophores per frame using emitter data
molecules_per_frame = zeros(Int, nframes)
for emitter in smld_noisy.emitters
    frame = Int(emitter.frame)
    if 1 <= frame <= nframes
        molecules_per_frame[frame] += 1
    end
end

# Plot molecules per frame
fig3 = Figure()
ax3 = Axis(fig3[1, 1], 
    title="Fluorophores per Frame", 
    xlabel="Frame", 
    ylabel="Number of Molecules",
    yticklabelsize=14)

# Plot the molecule count
lines!(ax3, 1:nframes, molecules_per_frame, linewidth=2, color=:blue)

# Add annotation for rates
text!(ax3, 0.05, 0.95, 
    text="k_off = $(round(k_off, digits=2)) Hz\nk_on = $(round(k_on, digits=2)) Hz\nTotal molecules: $(length(smld_true.emitters))",
    align=(:left, :top), 
    space=:relative, 
    fontsize=14)

display(fig3)

# Calculate statistics
avg_molecules = mean(molecules_per_frame)
max_molecules = maximum(molecules_per_frame)
total_activations = sum(molecules_per_frame)

println("\nMolecule activation statistics:")
println("- Average molecules per frame: $(round(avg_molecules, digits=2))")
println("- Maximum molecules per frame: $max_molecules")
println("- Total molecule activations: $total_activations")
