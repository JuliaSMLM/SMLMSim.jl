"""
Generate MP4 video showing sCMOS camera artifacts in diffusion simulation.

Parameters match user request:
- D = 0.1 μm²/s
- 0.1 μm pixels
- 100 fps
- 50 frames (reduced from 200 for speed)
- 1 molecule/μm²
- 200 photons per frame
"""

using SMLMSim
using SMLMData
using MicroscopePSFs
using Statistics
using VideoIO
using ColorTypes
using FixedPointNumbers

println("="^80)
println("sCMOS Video Generation")
println("="^80)

# Parameters
D = 0.1
pixel_size = 0.1
fps = 100
n_frames = 50  # Reduced for demo speed
density = 1.0
photons = 200.0
box_size = 5.0  # 5x5 μm for speed

println("\nParameters:")
println("  D = $D μm²/s")
println("  Pixel size = $pixel_size μm")
println("  FPS = $fps")
println("  Frames = $n_frames")
println("  Density = $density /μm²")
println("  Photons = $photons")
println("  Box size = $(box_size)x$(box_size) μm")

n_pixels = Int(ceil(box_size / pixel_size))
println("  Image size: $(n_pixels)x$(n_pixels) pixels")

# Create extreme sCMOS camera
println("\nCreating sCMOS camera with extreme artifacts...")
readnoise_map = zeros(n_pixels, n_pixels)
offset_map = zeros(n_pixels, n_pixels)
gain_map = ones(n_pixels, n_pixels)

for i in 1:n_pixels
    for j in 1:n_pixels
        # Stripe pattern in readnoise
        readnoise_map[i, j] = 3.0 + 7.0 * abs(sin(2π * i / 5))

        # Strong gradient + stripes in offset
        gradient = 80.0 * (i / n_pixels)
        stripe = xor(mod(div(i-1, 3), 2), mod(div(j-1, 3), 2)) == 1 ? 40.0 : 0.0
        offset_map[i, j] = 20.0 + gradient + stripe

        # Gain variation
        gain_map[i, j] = 1.0 + 0.4 * sin(2π * i / 8) * cos(2π * j / 8)
    end
end

camera_scmos = SCMOSCamera(
    collect(0:pixel_size:box_size),
    collect(0:pixel_size:box_size),
    offset_map,
    gain_map,
    readnoise_map,
    0.85
)

println("  Readnoise: $(round(minimum(readnoise_map), digits=1)) - $(round(maximum(readnoise_map), digits=1)) e⁻")
println("  Offset: $(round(minimum(offset_map), digits=1)) - $(round(maximum(offset_map), digits=1)) ADU")

# Run diffusion simulation
println("\nRunning diffusion simulation...")
params = DiffusionSMLMParams(
    density = density,
    box_size = box_size,
    dt = 0.001,
    t_max = n_frames / fps,
    camera_framerate = fps,
    camera_exposure = 1.0 / fps,
    diff_monomer = D
)

smld = simulate(params; camera=camera_scmos, photons=photons)
println("  Generated $(length(smld.emitters)) localizations")

# Generate images
println("\nGenerating images...")
psf = GaussianPSF(0.13)

println("  sCMOS images (with camera_noise=true)...")
images_scmos = gen_images(smld, psf, camera_noise=true, bg=5.0)

println("  Ideal images (for comparison)...")
camera_ideal = IdealCamera(n_pixels, n_pixels, pixel_size)
smld_ideal = BasicSMLD(smld.emitters, camera_ideal, smld.n_frames, smld.n_datasets)
images_ideal = gen_images(smld_ideal, psf, poisson_noise=true, bg=5.0)

# Create video
println("\nCreating video...")
mkpath("dev/outputs")

# Get intensity range for consistent scaling
vmin_s, vmax_s = quantile(vec(images_scmos), [0.001, 0.999])
vmin_i, vmax_i = quantile(vec(images_ideal), [0.001, 0.999])

println("  sCMOS range: $(round(vmin_s, digits=1)) - $(round(vmax_s, digits=1))")
println("  Ideal range: $(round(vmin_i, digits=1)) - $(round(vmax_i, digits=1))")

# Prepare frames (side-by-side comparison)
video_frames = []
for i in 1:min(n_frames, size(images_scmos, 3))
    frame_s = images_scmos[:, :, i]
    frame_i = images_ideal[:, :, i]

    # Normalize
    norm_s = clamp.((frame_s .- vmin_s) ./ (vmax_s - vmin_s), 0, 1)
    norm_i = clamp.((frame_i .- vmin_i) ./ (vmax_i - vmin_i), 0, 1)

    # Side by side
    combined = hcat(norm_s, norm_i)
    gray_frame = Gray{N0f8}.(combined)

    push!(video_frames, gray_frame)
end

# Encode MP4
output_mp4 = "dev/outputs/scmos_diffusion.mp4"
encoder_options = (color_range=2, crf=20, preset="medium")

println("  Encoding MP4 ($(length(video_frames)) frames at $fps fps)...")
open_video_out(output_mp4, video_frames[1], framerate=fps, encoder_options=encoder_options) do writer
    for (i, frame) in enumerate(video_frames)
        write(writer, frame)
        if mod(i, 10) == 0
            print(".")
        end
    end
end
println("\n  Saved: $output_mp4")

# Stats
println("\nVideo Statistics:")
println("  Format: $(n_pixels*2)x$(n_pixels) (sCMOS | Ideal)")
println("  Frames: $(length(video_frames))")
println("  Framerate: $fps fps")
println("  Duration: $(round(length(video_frames)/fps, digits=2)) seconds")

filesize_mb = stat(output_mp4).size / 1024^2
println("  File size: $(round(filesize_mb, digits=2)) MB")

println("\n" * "="^80)
println("✓ Video created successfully!")
println("  Open: $output_mp4")
println("  Left = sCMOS (with spatial artifacts)")
println("  Right = Ideal (Poisson noise only)")
println("="^80)
