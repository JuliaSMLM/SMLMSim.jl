# Generate camera images from SMLD and PSF

"""
    gen_images(smld::SMLD, psf::AbstractPSF; kwargs...) -> Array{T, 3} where T<:Real

Generate camera images from SMLD data using the specified PSF model.

# Arguments
- `smld::SMLD`: Single molecule localization data container
- `psf::AbstractPSF`: Point spread function model

# Keyword arguments
- `dataset::Int=1`: Dataset number to use from SMLD
- `frames=nothing`: Specific frames to generate (default: all frames in smld.n_frames)
- `frame_integration::Int=1`: Number of consecutive frames to integrate into one image
- `support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}}=Inf`: PSF support region size
- `sampling::Int=2`: Supersampling factor for PSF integration
- `threaded::Bool=true`: Enable multithreading for faster computation
- `bg::Float64=0.0`: Background signal level (photons per pixel)
- `poisson_noise::Bool=false`: Apply Poisson noise
- `camera_noise::Bool=false`: Apply camera read noise (Note: This feature is not yet implemented)

# Returns
- 3D array of camera images with dimensions [height, width, num_frames]
- The element type T matches the type of emitter.photons (typically Float64)
"""
function gen_images(smld::SMLD, psf::AbstractPSF; 
          dataset::Int=1,
          frames=nothing,
          frame_integration::Int=1, 
          support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}}=Inf,
          sampling::Int=2,
          threaded::Bool=true,
          bg::Float64=0.0,
          poisson_noise::Bool=false,
          camera_noise::Bool=false)
    
    # Filter for the specified dataset
    dataset_smld = @filter(smld, dataset == dataset)
    
    # Determine frames to process
    if frames === nothing
        # Use all frames from 1 to n_frames
        frames = 1:smld.n_frames
    end
    
    # Group frames for integration
    if frame_integration > 1
        frame_groups = [frames[i:min(i+frame_integration-1, length(frames))] 
                       for i in 1:frame_integration:length(frames)]
    else
        frame_groups = [[f] for f in frames]
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
    images = zeros(T, height, width, length(frame_groups))
    
    # Process each frame group
    for (i, frame_group) in enumerate(frame_groups)
        # Filter emitters for this frame group
        frame_smld = filter_frames(dataset_smld, frame_group)
        
        # Add background to all images
        images[:,:,i] .+= bg
        
        # Skip image generation if no emitters in this frame group
        if isempty(frame_smld.emitters)
            continue
        end
        
        # Generate image with all emitters at once
        img = integrate_pixels(
            psf, 
            smld.camera, 
            frame_smld.emitters; 
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
    
    # TODO: Apply camera noise if requested
    if camera_noise
        @warn "Camera noise model not yet implemented"
    end
    
    return images
end

"""
    gen_image(smld::SMLD, psf::AbstractPSF, frame::Int; kwargs...) -> Matrix{T} where T<:Real

Generate a single camera image for a specific frame from SMLD data.
See `gen_images` for full documentation of parameters.

# Returns
- A 2D camera image as Matrix{T} where T matches the type of emitter.photons
"""
function gen_image(smld::SMLD, psf::AbstractPSF, frame::Int; kwargs...)
    # Call gen_images with a single frame
    images = gen_images(smld, psf; frames=[frame], kwargs...)
    return images[:,:,1]
end
