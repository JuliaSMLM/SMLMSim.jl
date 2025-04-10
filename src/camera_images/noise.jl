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