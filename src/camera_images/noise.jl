# Camera noise models for SMLM imaging
# Distributions package is imported at the module level

"""
    poisson_noise(image::AbstractArray{T}) where T<:Real -> Array{Float64}

Apply Poisson noise to an image or image stack.

# Arguments
- `image::AbstractArray{T}`: Input image or image stack with values representing photon counts

# Returns
- Array with same dimensions as input, with Poisson noise applied to each pixel

# Details
This function creates a copy of the input array and applies Poisson noise 
to each pixel using the in-place poisson_noise! function.

Non-integer and negative values are handled specially:
- Non-integer values are accepted (treating them as expected photon counts)
- Negative values are clipped to zero before applying noise
- Zero values remain zero (as Poisson(0) always returns 0)

# Example
```julia
# Add Poisson noise to a clean image
clean_image = ones(100, 100) * 100.0  # 100 expected photons per pixel
noisy_image = poisson_noise(clean_image)
```
"""
function poisson_noise(image::AbstractArray{T}) where T<:Real
    # Create a copy of the input array with Float64 type
    result = similar(image, Float64)
    result .= image
    
    # Apply in-place function to the copy
    poisson_noise!(result)
    
    return result
end

"""
    poisson_noise!(image::AbstractArray{T}) where T<:Real -> nothing

Apply Poisson noise to an image or image stack in-place.

# Arguments
- `image::AbstractArray{T}`: Input image or image stack with values representing photon counts

# Returns
- `nothing`: The input array is modified in-place

# Details
Same as `poisson_noise`, but modifies the input array directly instead of creating a new one.
This can be more memory-efficient for large images or when processing multiple frames.

# Example
```julia
# Add Poisson noise to an image in-place
image = ones(100, 100) * 100.0  # 100 expected photons per pixel
poisson_noise!(image)  # image is modified in-place
```
"""
function poisson_noise!(image::AbstractArray{T}) where T<:Real
    # Apply Poisson noise to each element in-place
    for i in eachindex(image)
        # Clip negative values to 0
        λ = max(image[i], 0.0)

        # Sample from Poisson distribution
        if λ > 0.0
            image[i] = Float64(rand(Poisson(λ)))
        else
            image[i] = 0.0
        end
    end

    return nothing
end

"""
    scmos_noise(image::AbstractMatrix{T}, camera::SCMOSCamera) where T<:Real -> Matrix{Float64}

Apply realistic sCMOS camera noise model to an image.

# Arguments
- `image::AbstractMatrix{T}`: Input image with values representing photon counts
- `camera::SCMOSCamera`: sCMOS camera with calibration parameters (offset, gain, readnoise, qe)

# Returns
- Matrix with same dimensions as input, with full sCMOS noise model applied

# Details
The sCMOS noise model applies the following transformations in order:
1. Quantum efficiency: Convert photons to photoelectrons
2. Poisson noise: Shot noise on photoelectrons
3. Read noise: Gaussian noise per pixel
4. Gain: Convert electrons to ADU (analog-to-digital units)
5. Offset: Add dark level

The process simulates the physical detection chain:
- QE: Not all photons generate photoelectrons
- Poisson: Fundamental shot noise
- Read noise: Electronic noise from amplifier/ADC
- Gain & Offset: Conversion to digital units

# Example
```julia
# Create an sCMOS camera with 1.6 e⁻ read noise
camera = SCMOSCamera(128, 128, 0.1, 1.6)

# Apply realistic noise to a clean image
clean_image = ones(128, 128) * 100.0  # 100 photons per pixel
noisy_image = scmos_noise(clean_image, camera)
```
"""
function scmos_noise(image::AbstractMatrix{T}, camera::SCMOSCamera) where T<:Real
    height, width = size(image)
    result = similar(image, Float64)

    # Apply sCMOS noise model to each pixel
    for j in 1:width
        for i in 1:height
            # Get pixel-specific calibration parameters
            qe = SMLMData.get_qe(camera, i, j)
            gain = SMLMData.get_gain(camera, i, j)
            offset = SMLMData.get_offset(camera, i, j)
            readnoise = SMLMData.get_readnoise(camera, i, j)

            # 1. Apply quantum efficiency: photons -> photoelectrons
            photons = max(image[i, j], 0.0)
            photoelectrons_mean = photons * qe

            # 2. Apply Poisson noise (shot noise)
            if photoelectrons_mean > 0.0
                photoelectrons = Float64(rand(Poisson(photoelectrons_mean)))
            else
                photoelectrons = 0.0
            end

            # 3. Add read noise (Gaussian, in electrons)
            if readnoise > 0.0
                photoelectrons += randn() * readnoise
            end

            # 4. Convert to ADU using gain: electrons -> ADU
            adu = photoelectrons / gain

            # 5. Add offset (dark level)
            result[i, j] = adu + offset
        end
    end

    return result
end

"""
    scmos_noise!(image::AbstractMatrix{T}, camera::SCMOSCamera) where T<:Real -> nothing

Apply realistic sCMOS camera noise model to an image in-place.

# Arguments
- `image::AbstractMatrix{T}`: Input image with values representing photon counts (modified in-place)
- `camera::SCMOSCamera`: sCMOS camera with calibration parameters

# Returns
- `nothing`: The input array is modified in-place

# Details
Same as `scmos_noise`, but modifies the input array directly instead of creating a new one.

# Example
```julia
camera = SCMOSCamera(128, 128, 0.1, 1.6)
image = ones(128, 128) * 100.0
scmos_noise!(image, camera)  # image is modified in-place
```
"""
function scmos_noise!(image::AbstractMatrix{T}, camera::SCMOSCamera) where T<:Real
    height, width = size(image)

    # Apply sCMOS noise model to each pixel in-place
    for j in 1:width
        for i in 1:height
            # Get pixel-specific calibration parameters
            qe = SMLMData.get_qe(camera, i, j)
            gain = SMLMData.get_gain(camera, i, j)
            offset = SMLMData.get_offset(camera, i, j)
            readnoise = SMLMData.get_readnoise(camera, i, j)

            # 1. Apply quantum efficiency: photons -> photoelectrons
            photons = max(image[i, j], 0.0)
            photoelectrons_mean = photons * qe

            # 2. Apply Poisson noise (shot noise)
            if photoelectrons_mean > 0.0
                photoelectrons = Float64(rand(Poisson(photoelectrons_mean)))
            else
                photoelectrons = 0.0
            end

            # 3. Add read noise (Gaussian, in electrons)
            if readnoise > 0.0
                photoelectrons += randn() * readnoise
            end

            # 4. Convert to ADU using gain: electrons -> ADU
            adu = photoelectrons / gain

            # 5. Add offset (dark level)
            image[i, j] = adu + offset
        end
    end

    return nothing
end