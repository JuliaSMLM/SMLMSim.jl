# test_stack_viewer.jl - Test the image stack viewer with synthetic data
using Pkg
Pkg.activate("dev")
using Revise

using GLMakie
include("stack_viewer.jl")

"""
Generate a 3D test image stack with various patterns.

# Arguments
- `height`: Image height in pixels
- `width`: Image width in pixels
- `n_frames`: Number of frames in the stack
- `noise_level`: Amount of background noise to add (0-1)

# Returns
- 3D array representing the image stack (height × width × frames)
"""
function generate_test_stack(height=256, width=256, n_frames=50, noise_level=0.1)
    # Create empty stack
    stack = zeros(Float32, height, width, n_frames)

    # Add background noise
    stack .+= noise_level * rand(Float32, height, width, n_frames)

    # Add moving Gaussian spots
    for frame in 1:n_frames
        # Calculate position based on frame (circular motion)
        angle = 2π * frame / n_frames
        center_x = width / 2 + (width / 3) * cos(angle)
        center_y = height / 2 + (height / 3) * sin(angle)

        # Create Gaussian spot
        for y in 1:height, x in 1:width
            # Distance from center
            dist = sqrt((x - center_x)^2 + (y - center_y)^2)
            # Gaussian intensity with sigma based on frame number
            sigma = 5 + 10 * (frame / n_frames)
            intensity = exp(-dist^2 / (2 * sigma^2))
            stack[y, x, frame] += intensity
        end

        # Add some static spots for reference
        for spot in 1:5
            spot_x = width * (0.2 + 0.6 * (spot / 5))
            spot_y = height * 0.8

            for y in 1:height, x in 1:width
                dist = sqrt((x - spot_x)^2 + (y - spot_y)^2)
                intensity = 0.8 * exp(-dist^2 / (2 * 6^2))
                stack[y, x, frame] += intensity
            end
        end
    end

    return stack
end

println("Generating test image stack...")
img_stack = generate_test_stack(256, 256, 50, 0.05)

println("Stack dimensions: $(size(img_stack))")
println("Min value: $(minimum(img_stack)), Max value: $(maximum(img_stack))")

println("Launching stack viewer...")
fig = view_stack(img_stack, title="Test Stack Viewer")

# The viewer is now interactive
println("Viewer is now interactive. Close the window to exit.")

display(fig)
