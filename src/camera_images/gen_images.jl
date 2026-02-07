# Generate camera images from SMLD and PSF

"""
    gen_images(smld::SMLD, psf::AbstractPSF; kwargs...) -> Tuple{Array{T, 3}, ImageInfo} where T<:Real

Generate camera images from SMLD data using the specified PSF model.

# Arguments
- `smld::SMLD`: Single molecule localization data container
- `psf::AbstractPSF`: Point spread function model

# Keyword arguments
- `dataset::Int=1`: Dataset number to use from SMLD
- `frames=nothing`: Specific frames to generate (default: all frames in smld.n_frames)
- `support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}}=Inf`: PSF support region size:
  - `Inf` (default): Calculate PSF over the entire image (most accurate but slowest)
  - `Real`: Circular region with specified radius (in microns) around each emitter
  - `Tuple{<:Real,<:Real,<:Real,<:Real}`: Explicit region as (xmin, xmax, ymin, ymax) in microns
- `sampling::Int=2`: Supersampling factor for PSF integration
- `threaded::Bool=true`: Enable multithreading for faster computation
- `bg::Float64=0.0`: Background signal level (photons per pixel)
- `poisson_noise::Bool=false`: Apply Poisson noise only (for simple shot noise)
- `camera_noise::Bool=false`: Apply full camera noise model (requires SCMOSCamera)
  - For SCMOSCamera: applies QE, Poisson, read noise, gain, and offset
  - For IdealCamera: ignored (use poisson_noise instead)

# Returns
- `Tuple{Array{T,3}, ImageInfo}`: (images, info)
    - images: 3D array of camera images with dimensions [height, width, num_frames]
    - info: ImageInfo containing timing and image statistics

# Performance Note
For the `support` parameter, using a finite radius (typically 3-5× the PSF width)
provides a good balance between accuracy and performance. For example, with a PSF
width of 0.15μm, a support radius of 0.5-1.0μm is usually sufficient.
"""
function gen_images(smld::SMLD, psf::AbstractPSF;
          dataset::Int=1,
          frames=nothing,
          support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}}=Inf,
          sampling::Int=2,
          threaded::Bool=true,
          bg::Float64=0.0,
          poisson_noise::Bool=false,
          camera_noise::Bool=false)

    start_time = time_ns()

    # Filter for the specified dataset
    dataset_smld = @filter(smld, dataset == dataset)

    # Determine frames to process
    if frames === nothing
        # Use all frames from 1 to n_frames
        frames = 1:smld.n_frames
    end

    # Get camera dimensions from SMLD
    camera = smld.camera
    width = length(camera.pixel_edges_x) - 1
    height = length(camera.pixel_edges_y) - 1

    # Determine the element type from emitters (if any)
    if isempty(dataset_smld.emitters)
        T = Float64 # default if no emitters
    else
        T = typeof(dataset_smld.emitters[1].photons)
    end

    # Pre-allocate output array
    images = zeros(T, height, width, length(frames))

    # Track total photons
    n_photons_total = 0.0

    # Process each frame individually
    for (i, frame_num) in enumerate(frames)
        # Filter emitters for this frame
        frame_emitters = filter(e -> e.frame == frame_num, dataset_smld.emitters)

        # Add background to all images
        images[:,:,i] .+= bg

        # Skip image generation if no emitters in this frame
        if isempty(frame_emitters)
            continue
        end

        # Accumulate photon count
        n_photons_total += sum(e.photons for e in frame_emitters)

        # Generate image with all emitters at once
        img = integrate_pixels(
            psf,
            smld.camera,
            frame_emitters;
            support=support,
            sampling=sampling,
            threaded=threaded
        )

        # Add image to background
        images[:,:,i] .+= img
    end

    # Apply Poisson noise if requested
    if poisson_noise
        # Loop through each pixel and apply Poisson noise
        for i in eachindex(images)
            # Apply Poisson noise directly
            λ = max(images[i], 0.0)
            if λ > 0.0
                images[i] = Float64(rand(Poisson(λ)))
            else
                images[i] = 0.0
            end
        end
    end

    # Apply camera noise if requested
    if camera_noise
        if camera isa SCMOSCamera
            # Apply full sCMOS noise model (QE, Poisson, read noise, gain, offset) to each frame
            for frame_idx in 1:size(images, 3)
                frame = @view images[:, :, frame_idx]
                scmos_noise!(frame, camera)
            end
        else
            # For IdealCamera, camera_noise flag is ignored (IdealCamera is Poisson-only)
            @warn "Camera noise only supported for SCMOSCamera. Use poisson_noise=true for IdealCamera." maxlog=1
        end
    end

    elapsed_s = (time_ns() - start_time) / 1e9

    # Build ImageInfo
    info = ImageInfo(
        elapsed_s=elapsed_s,
        backend=:cpu,
        device_id=-1,
        frames_generated=length(frames),
        n_photons_total=n_photons_total,
        output_size=(height, width, length(frames))
    )

    return images, info
end

"""
    gen_image(smld::SMLD, psf::AbstractPSF, frame::Int; kwargs...) -> Tuple{Matrix{T}, ImageInfo} where T<:Real

Generate a single camera image for a specific frame from SMLD data.
See `gen_images` for full documentation of parameters.

# Returns
- `Tuple{Matrix{T}, ImageInfo}`: (image, info)
    - image: 2D camera image as Matrix{T}
    - info: ImageInfo containing timing and image statistics
"""
function gen_image(smld::SMLD, psf::AbstractPSF, frame::Int; kwargs...)
    # Call gen_images with a single frame
    images, info = gen_images(smld, psf; frames=[frame], kwargs...)
    return images[:,:,1], info
end
