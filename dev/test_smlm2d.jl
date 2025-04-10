# test_smlm2d.jl - Generate 2D SMLM data with Nmer pattern and view using stack viewer

using Pkg
Pkg.activate("dev")
using Revise

using SMLMSim
using MicroscopePSFs
using GLMakie
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
    nframes = 100,     # Generate 100 frames
    framerate = 20.0,  # 20 frames per second
    ndims = 2,         # 2D simulation
    zrange = [-1.0, 1.0]  # Not used for 2D, but required parameter
)

println("Running SMLM simulation with hexamer pattern...")

# Create a fluorophore model with custom blinking kinetics

# Number of blinking events per fluorophore
n_blinks_per_fluor = 5 
total_time = params.nframes / params.framerate  # Total time in seconds

k_on = 1/(total_time / n_blinks_per_fluor)  # On state rate (Hz)
k_off = params.framerate  # Off state rate (Hz)
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
    bg=10.0,             # Background photons
    poisson_noise=true   # Add realistic Poisson noise
)

println("Image stack generated with size: $(size(images))")

# Launch the stack viewer
println("Launching stack viewer...")
fig2 = view_stack(images, title="SMLM 2D Simulation (Hexamer Pattern)")

println("Stack viewer is now interactive. Close the window to exit.")

# # Check fluorophores per frame 

# fluor_per_frame = vec(sum(images, dims = [1,2]))
# fig3 = Figure()
# ax3 = Axis(fig3[1, 1], title="Fluorophores per Frame", xlabel="Frame", ylabel="Fluorophores")
# lines!(ax3, 1:length(fluor_per_frame), fluor_per_frame)
# display(fig3)
