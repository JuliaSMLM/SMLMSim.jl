"""
Demonstration of sCMOS camera noise artifacts in diffusion simulation.

Creates an extreme sCMOS camera with spatially-varying noise characteristics
to make artifacts clearly visible in the output video.

Parameters:
- D = 0.1 μm²/s diffusion coefficient
- 0.1 μm pixel size
- 100 fps (10 ms exposure)
- 200 frames (2 seconds total)
- 1 molecule/μm² density
- 200 photons per localization
"""

using SMLMSim
using SMLMData
using Printf

println("=" ^ 80)
println("sCMOS Camera Noise Demonstration")
println("=" ^ 80)

# Simulation parameters
D = 0.1          # μm²/s
pixel_size = 0.1 # μm
fps = 100        # frames per second
n_frames = 50    # total frames (reduced for faster demo)
density = 1.0    # molecules per μm²
photons = 200.0  # photons per localization
box_size = 10.0  # μm (10x10 μm field of view)

println("\nSimulation Parameters:")
println("  Diffusion coefficient: $D μm²/s")
println("  Pixel size: $pixel_size μm")
println("  Frame rate: $fps fps")
println("  Total frames: $n_frames")
println("  Density: $density molecules/μm²")
println("  Photons: $photons per localization")
println("  Field of view: $(box_size)x$(box_size) μm²")

# Calculate camera dimensions
n_pixels = Int(ceil(box_size / pixel_size))
println("\nCamera: $(n_pixels)x$(n_pixels) pixels")

# Create extreme sCMOS camera with spatially varying noise
println("\nCreating sCMOS camera with extreme spatially-varying noise...")

# Create base camera
camera_exposure = 1.0 / fps  # seconds
readnoise_base = 5.0  # e⁻ rms

# Create per-pixel readnoise map with dramatic variation
# Make it have stripes and hot/cold regions
readnoise_map = zeros(n_pixels, n_pixels)
for i in 1:n_pixels
    for j in 1:n_pixels
        # Create stripe pattern in readnoise
        stripe = 1.0 + 0.8 * sin(2π * i / 10)

        # Add hot corners
        dist_corner = min(
            sqrt((i-1)^2 + (j-1)^2),
            sqrt((i-n_pixels)^2 + (j-1)^2),
            sqrt((i-1)^2 + (j-n_pixels)^2),
            sqrt((i-n_pixels)^2 + (j-n_pixels)^2)
        )
        hot_factor = 1.0 + 3.0 * exp(-dist_corner / 5)

        readnoise_map[i, j] = readnoise_base * stripe * hot_factor
    end
end

# Create per-pixel offset map with dramatic gradients
offset_map = zeros(n_pixels, n_pixels)
for i in 1:n_pixels
    for j in 1:n_pixels
        # Strong gradient across the sensor
        gradient = 50.0 + 30.0 * (i / n_pixels) + 20.0 * (j / n_pixels)

        # Add checkerboard pattern
        checker = xor(mod(div(i-1, 5), 2), mod(div(j-1, 5), 2)) == 1 ? 10.0 : -10.0

        offset_map[i, j] = gradient + checker
    end
end

# Create gain map with variation
gain_map = zeros(n_pixels, n_pixels)
for i in 1:n_pixels
    for j in 1:n_pixels
        # Gain varies from 0.8 to 1.2
        gain_map[i, j] = 1.0 + 0.2 * sin(2π * i / 15) * cos(2π * j / 15)
    end
end

# Create sCMOS camera with these maps
camera_scmos = SCMOSCamera(
    collect(0:pixel_size:box_size),  # pixel_edges_x
    collect(0:pixel_size:box_size),  # pixel_edges_y
    offset_map,
    gain_map,
    readnoise_map,
    0.95  # constant QE (could also make this spatially varying)
)

println("  Read noise range: $(minimum(readnoise_map)) - $(maximum(readnoise_map)) e⁻")
println("  Offset range: $(minimum(offset_map)) - $(maximum(offset_map)) ADU")
println("  Gain range: $(minimum(gain_map)) - $(maximum(gain_map))")

# Set up diffusion simulation parameters
println("\nRunning diffusion simulation...")
params = DiffusionSMLMParams(
    density = density,
    box_size = box_size,
    dt = 0.001,  # 1 ms timesteps for simulation
    t_max = n_frames / fps,  # Total time in seconds
    camera_framerate = fps,
    camera_exposure = camera_exposure,
    diff_monomer = D
)

# Run simulation
smld = simulate(params; camera=camera_scmos, photons=photons)

println("  Generated $(length(smld.emitters)) localizations")
println("  Frames: $(smld.n_frames)")

# Generate images with sCMOS camera noise
println("\nGenerating images with sCMOS camera noise...")
using MicroscopePSFs
psf = GaussianPSF(0.13)  # 130 nm PSF width

images_scmos = gen_images(smld, psf, camera_noise=true, bg=10.0)
println("  Image stack size: $(size(images_scmos))")

# Also generate ideal images for comparison
println("\nGenerating ideal comparison images...")
camera_ideal = IdealCamera(n_pixels, n_pixels, pixel_size)
smld_ideal = BasicSMLD(smld.emitters, camera_ideal, smld.n_frames, smld.n_datasets)
images_ideal = gen_images(smld_ideal, psf, poisson_noise=true, bg=10.0)

# Calculate statistics
println("\nImage Statistics:")
println("  sCMOS - Mean: $(round(mean(images_scmos), digits=2)), Std: $(round(std(images_scmos), digits=2))")
println("  sCMOS - Range: $(round(minimum(images_scmos), digits=2)) - $(round(maximum(images_scmos), digits=2))")
println("  Ideal - Mean: $(round(mean(images_ideal), digits=2)), Std: $(round(std(images_ideal), digits=2))")
println("  Ideal - Range: $(round(minimum(images_ideal), digits=2)) - $(round(maximum(images_ideal), digits=2))")

# Save results
println("\nSaving results to dev/outputs/...")

# Save as HDF5 for analysis
using HDF5
output_file = "dev/outputs/scmos_demo.h5"
h5open(output_file, "w") do file
    write(file, "images_scmos", images_scmos)
    write(file, "images_ideal", images_ideal)
    write(file, "readnoise_map", readnoise_map)
    write(file, "offset_map", offset_map)
    write(file, "gain_map", gain_map)

    # Save parameters as attributes
    attrs(file)["D"] = D
    attrs(file)["pixel_size"] = pixel_size
    attrs(file)["fps"] = fps
    attrs(file)["n_frames"] = n_frames
    attrs(file)["density"] = density
    attrs(file)["photons"] = photons
end
println("  Saved: $output_file")

# Create visualization and save as mp4
println("\nCreating video visualization...")

# We'll need to check if we have video encoding capability
# Try using VideoIO and FFMPEG
try
    using VideoIO
    using ColorTypes
    using FixedPointNumbers

    output_video = "dev/outputs/scmos_demo.mp4"

    # Prepare frames for video (normalize and create side-by-side comparison)
    # We'll show: sCMOS | Ideal | Difference
    video_frames = []

    for i in 1:n_frames
        frame_scmos = images_scmos[:, :, i]
        frame_ideal = images_ideal[:, :, i]

        # Normalize each to 0-1 for display
        vmin_s, vmax_s = quantile(vec(images_scmos), [0.01, 0.99])
        vmin_i, vmax_i = quantile(vec(images_ideal), [0.01, 0.99])

        norm_scmos = clamp.((frame_scmos .- vmin_s) ./ (vmax_s - vmin_s), 0, 1)
        norm_ideal = clamp.((frame_ideal .- vmin_i) ./ (vmax_i - vmin_i), 0, 1)

        # Difference (normalized separately)
        diff = frame_scmos - frame_ideal
        vmin_d, vmax_d = quantile(vec(diff), [0.01, 0.99])
        norm_diff = clamp.((diff .- vmin_d) ./ (vmax_d - vmin_d), 0, 1)

        # Create side-by-side (convert to Gray)
        combined = hcat(norm_scmos, norm_ideal, norm_diff)
        gray_frame = Gray{N0f8}.(combined)

        push!(video_frames, gray_frame)
    end

    # Encode video
    encoder_options = (color_range=2, crf=23, preset="medium")
    framerate = fps
    open_video_out(output_video, video_frames[1], framerate=framerate, encoder_options=encoder_options) do writer
        for frame in video_frames
            write(writer, frame)
        end
    end

    println("  Saved video: $output_video")
    println("  Video shows: sCMOS | Ideal | Difference")

catch e
    println("  Warning: Could not create video (VideoIO not available)")
    println("  Error: $e")
    println("  Data saved to HDF5 file for manual visualization")
end

println("\n" * "=" ^ 80)
println("Demo complete! Check dev/outputs/ for results.")
println("=" ^ 80)
