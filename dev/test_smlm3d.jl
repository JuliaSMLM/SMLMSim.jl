# test_smlm3d.jl - Generate 3D SMLM data with Nmer pattern and view using stack viewer

using Pkg
Pkg.activate("dev")
using Revise

using SMLMSim
using MicroscopePSFs
using GLMakie
using Statistics
include("stack_viewer.jl")

println("Setting up 3D simulation parameters...")

# Create a camera with 100nm pixels
pixelsize = 0.1  # 100nm pixels
camera = IdealCamera(1:256, 1:256, pixelsize)

# Create a 3D pattern - 8 molecules arranged in a circle with 300nm diameter
pattern = Nmer3D(n=8, d=0.3)  # 8 molecules in a 300nm circle at z=0

# Create static SMLM parameters for 3D simulation
params = StaticSMLMParams(
    ρ = 0.5,           # 0.5 patterns per square micron (less dense than 2D example)
    σ_psf = 0.13,      # 130nm PSF width (realistic for visible light)
    minphotons = 100,  # Minimum photon count for detection
    ndatasets = 1,     # Generate 1 dataset
    nframes = 100,     # Generate 100 frames
    framerate = 10.0,  # 10 frames per second
    ndims = 3,         # 3D simulation
    zrange = [-2.0, 2.0]  # 4μm axial range
)

println("Running 3D SMLM simulation with Nmer3D pattern...")

# Create a fluorophore model with custom blinking kinetics
nframes = params.nframes
framerate = params.framerate
total_time = nframes / framerate  # Total time in seconds

# Kinetics parameters
k_off = framerate        # Off state rate (Hz) - corresponds to 1 frame duration
n_blinks_per_fluor = 3   # Desired number of events per 100 frames
k_on = n_blinks_per_fluor / total_time  # On state rate (Hz)
photons_per_frame = 2000  # More photons for 3D imaging
photons = photons_per_frame * framerate  # Photon emission rate (Hz)

println("Fluorophore parameters:")
println("- Off rate: $k_off Hz (average on time: $(1/k_off) seconds)")
println("- On rate: $k_on Hz (average off time: $(1/k_on) seconds)")
println("- Expected events per $nframes frames: $n_blinks_per_fluor")
println("- Photons per frame: $photons_per_frame")

fluor = GenericFluor(
    photons=photons,  # Photon emission rate (Hz)
    k_off=k_off,      # Off rate (on→off) in Hz
    k_on=k_on         # On rate (off→on) in Hz
)

# Run simulation with custom pattern and fluorophore model
smld_true, smld_model, smld_noisy = simulate(
    params;
    pattern=pattern,
    molecule=fluor,
    camera=camera
)

println("3D SMLM simulation complete")
println("- True emitters: $(length(smld_true.emitters))")
println("- Model emitters after blinking: $(length(smld_model.emitters))")
println("- Noisy emitters with localization uncertainty: $(length(smld_noisy.emitters))")

# Extract coordinates
x_true = [e.x for e in smld_true.emitters]
y_true = [e.y for e in smld_true.emitters]
z_true = [e.z for e in smld_true.emitters]

x_noisy = [e.x for e in smld_noisy.emitters]
y_noisy = [e.y for e in smld_noisy.emitters]
z_noisy = [e.z for e in smld_noisy.emitters]
frames = [e.frame for e in smld_noisy.emitters]
photon_counts = [e.photons for e in smld_noisy.emitters]

# Plot 3D visualization
fig1 = Figure(size=(1200, 600))

# 3D scatter plots for true and noisy positions using LScene with scenekw
lscene1 = LScene(fig1[1, 1], scenekw=(show_axis=true,))
scatter!(lscene1, x_true, y_true, z_true, color=:blue, markersize=8)

# Add manual axis labels for 3D scene
Label(fig1[1, 1, Bottom()], "X (μm)", padding=(0, 0, 0, 5))
Label(fig1[1, 1, Left()], "Y (μm)", rotation=π/2, padding=(0, 5, 0, 0))

# Add title above the scene
Label(fig1[1, 1, Top()], "True Emitter Positions", fontsize=16)

lscene2 = LScene(fig1[1, 2], scenekw=(show_axis=true,))
scatter!(lscene2, x_noisy, y_noisy, z_noisy, color=z_noisy, colormap=:viridis, 
         markersize=8)

# Add manual axis labels for the second 3D scene
Label(fig1[1, 2, Bottom()], "X (μm)", padding=(0, 0, 0, 5))
Label(fig1[1, 2, Left()], "Y (μm)", rotation=π/2, padding=(0, 5, 0, 0))

# Add title above the scene
Label(fig1[1, 2, Top()], "Noisy Emitter Positions", fontsize=16)

Colorbar(fig1[1, 3], colormap=:viridis, label="Z position (μm)")

# Create 2D projections
ax3 = Axis(fig1[2, 1], 
    title="XY Projection (True)",
    xlabel="X (μm)", ylabel="Y (μm)",
    aspect=DataAspect())
scatter!(ax3, x_true, y_true, color=z_true, markersize=8)

ax4 = Axis(fig1[2, 2], 
    title="XY Projection (Noisy)",
    xlabel="X (μm)", ylabel="Y (μm)",
    aspect=DataAspect())
scatter!(ax4, x_noisy, y_noisy, color=z_noisy, colormap=:viridis, markersize=8)

# Link XY planes
linkaxes!(ax3, ax4)

# Use the simplest approach to link cameras - use the same camera for both scenes
lscene2.scene.camera = lscene1.scene.camera

# Function to reset all views (3D scenes and 2D plots)
function reset_all_views!()
    # Reset 2D axes
    autolimits!(ax3)
    autolimits!(ax4)
    
    # Reset 3D scene view using built-in zoom! function
    zoom!(lscene1.scene, 1.0)
end

# Add info banner at the top
Label(fig1[0, 1:2], "3D SMLM Simulation Results", 
      fontsize=16, font=:bold, padding=(0, 0, 5, 0))

# Add a reset button at the bottom with more prominence
fig1[3, 1:2] = GridLayout()
reset_btn = Button(fig1[3, 1:2], label="Reset Camera View", 
                  buttoncolor=:lightblue, fontsize=14)
on(reset_btn.clicks) do n
    reset_all_views!()
end

# Add rotation controls to explain to the user they can rotate both scenes
Label(fig1[0, 1:2], "Rotate either 3D scene to view from different angles (both scenes rotate together)", 
      fontsize=12, tellwidth=false, padding=(0, 0, 10, 0))

display(fig1)

# Create a PSF model for generating camera images
println("Creating 3D PSF model...")

# For 3D imaging, use astigmatic PSF or other 3D PSF model
zc = ZernikeCoefficients(15)
zc.phase[6] = 0.5 # Zernike coefficient for astigmatism
psf_raw = ScalarPSF(1.2, 0.6, 1.33; zernike_coeffs = zc)  # 130nm PSF width
x_range = y_range = range(-2.0, 2.0, step=0.05)  # 41×41 lateral grid, 4μm range
z_range = range(-1.0, 1.0, step=0.1)            # 21 z-planes, 2μm range
psf = SplinePSF(psf_raw, x_range, y_range, z_range)

# Generate camera images
println("Generating camera images...")
images = gen_images(smld_noisy, psf;
    support=1.0,        # PSF support radius in μm
    bg=2.0,             # Background photons
    poisson_noise=true  # Add realistic Poisson noise
)

println("Image stack generated with size: $(size(images))")

# Launch the stack viewer
println("Launching stack viewer...")
fig2 = view_stack(images, title="SMLM 3D Simulation (Nmer Pattern)")

# Run the stack viewer
println("Stack viewer is now interactive. Close the window to exit.")

println("\nSimulation and analysis complete.")