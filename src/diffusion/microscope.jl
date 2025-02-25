"""
    gen_image(psf::AbstractPSF, system::DiffusingMoleculeSystem, frame::Int;
              photons::Float64=1000.0, bg::Float64=5.0, poisson_noise::Bool=true)

Generate a microscope image from a diffusing molecule system.

# Arguments
- `psf::AbstractPSF`: Point spread function model (e.g., Airy2D, Gaussian2D)
- `system::DiffusingMoleculeSystem`: The molecular system state
- `frame::Int`: Frame number to generate

# Keyword Arguments
- `photons::Float64=1000.0`: Base photon count for emitters
- `bg::Float64=5.0`: Background photon count per pixel
- `poisson_noise::Bool=true`: Whether to add Poisson noise to final image

# Returns
- `Array{Float64,2}`: Generated image with dimensions matching system camera

# Example
```julia
# Generate an image from a system
psf = Gaussian2D(0.15)  # 150nm PSF width
image = gen_image(psf, system, 1; 
                 photons=2000.0, 
                 bg=10.0,
                 poisson_noise=true)
```
"""
function gen_image(psf::AbstractPSF, system::DiffusingMoleculeSystem, frame::Int;
                  photons::Float64=1000.0, bg::Float64=5.0, poisson_noise::Bool=true)
    
    # Input validation
    if frame < 1
        throw(ArgumentError("Frame number must be positive"))
    end
    
    if photons <= 0
        throw(ArgumentError("Photon count must be positive"))
    end
    
    if bg < 0
        throw(ArgumentError("Background level must be non-negative"))
    end
    
    # Get camera dimensions from pixel edges
    ny = length(system.camera.pixel_edges_y) - 1
    nx = length(system.camera.pixel_edges_x) - 1
    image = zeros(Float64, ny, nx)
    
    # Process each molecule
    for mol in system.molecules
        # Set photon count for this molecule
        mol.emitter.photons = photons
        
        # Skip if molecule isn't emitting (e.g. part of a non-fluorescent complex)
        mol.photons <= 0 && continue
        
        # Integrate PSF over pixels for this emitter
        pixel_intensities = integrate_pixels(psf, system.camera, mol.emitter)
        
        # Add to image
        image .+= pixel_intensities
    end
    
    # Add background
    image .+= bg
    
    # Apply Poisson noise if requested
    if poisson_noise
        for idx in eachindex(image)
            image[idx] = rand(Poisson(image[idx]))
        end
    end
    
    return image
end

"""
    gen_image_sequence(psf::AbstractPSF, systems::Vector{<:DiffusingMoleculeSystem};
                      photons::Float64=1000.0, bg::Float64=5.0,
                      frame_integration::Int=1, poisson_noise::Bool=true)

Generate a sequence of microscope images from diffusing molecule system states.

# Arguments
- `psf::AbstractPSF`: Point spread function model
- `systems::Vector{<:DiffusingMoleculeSystem}`: Sequence of system states

# Keyword Arguments
- `photons::Float64=1000.0`: Base photon count for emitters
- `bg::Float64=5.0`: Background photon count per pixel
- `frame_integration::Int=1`: Number of simulation frames to integrate per output frame
- `poisson_noise::Bool=true`: Whether to add Poisson noise

# Returns
- `Array{Float64,3}`: Stack of generated images [ny, nx, frames]

# Example
```julia
# Simulate diffusion
params = SmoluchowskiParams()
systems = simulate(params)

# Generate images with frame integration
psf = Gaussian2D(0.15)
image_stack = gen_image_sequence(psf, systems;
                                photons=1000.0,
                                bg=5.0,
                                frame_integration=10,
                                poisson_noise=true)
```
"""
function gen_image_sequence(psf::AbstractPSF, systems::Vector{<:DiffusingMoleculeSystem};
                          photons::Float64=1000.0, bg::Float64=5.0,
                          frame_integration::Int=1, poisson_noise::Bool=true)
    
    # Input validation
    if isempty(systems)
        throw(ArgumentError("Systems vector cannot be empty"))
    end
    
    if photons <= 0
        throw(ArgumentError("Photon count must be positive"))
    end
    
    if bg < 0
        throw(ArgumentError("Background level must be non-negative"))
    end
    
    if frame_integration < 1
        throw(ArgumentError("Frame integration must be at least 1"))
    end
    
    # Get dimensions from first system
    ny = length(systems[1].camera.pixel_edges_y) - 1
    nx = length(systems[1].camera.pixel_edges_x) - 1
    n_frames = ceil(Int, length(systems) / frame_integration)
    
    # Initialize image stack
    image_stack = zeros(Float64, ny, nx, n_frames)
    
    # Generate each frame, potentially integrating multiple states
    Threads.@threads for out_frame in 1:n_frames
        # Determine integration range
        start_idx = (out_frame - 1) * frame_integration + 1
        end_idx = min(start_idx + frame_integration - 1, length(systems))
        
        # Initialize integrated frame
        frame = zeros(Float64, ny, nx)
        
        # Integrate multiple system states
        for idx in start_idx:end_idx
            system = systems[idx]
            
            # Scale photons by integration time
            photon_scale = 1.0 / frame_integration
            for mol in system.molecules
                mol.emitter.photons = photons * photon_scale
            end
            
            # Generate and add frame without noise
            frame .+= gen_image(psf, system, idx, photons=photons, bg=0.0, poisson_noise=false)
        end
        
        # Add background after integration
        frame .+= bg
        
        # Apply noise to integrated frame if requested
        if poisson_noise
            for idx in eachindex(frame)
                frame[idx] = rand(Poisson(frame[idx]))
            end
        end
        
        image_stack[:, :, out_frame] = frame
    end
    
    return image_stack
end

"""
    simulate_and_image(params::SmoluchowskiParams, psf::AbstractPSF;
                      photons::Float64=1000.0, bg::Float64=5.0,
                      frame_integration::Int=1, poisson_noise::Bool=true)

Run a complete diffusion simulation and generate microscope images.

# Arguments
- `params::SmoluchowskiParams`: Simulation parameters
- `psf::AbstractPSF`: Point spread function model

# Keyword Arguments
- `photons::Float64=1000.0`: Base photon count for emitters
- `bg::Float64=5.0`: Background photon count per pixel
- `frame_integration::Int=1`: Number of simulation frames to integrate per output frame
- `poisson_noise::Bool=true`: Whether to add Poisson noise

# Returns
- `Tuple{Array{Float64,3}, Vector{<:DiffusingMoleculeSystem}}`: (images, systems)
  - `images`: Stack of generated images [ny, nx, frames]
  - `systems`: System states at each timestep

# Example
```julia
# Set up parameters and PSF
params = SmoluchowskiParams()
psf = Gaussian2D(0.15)

# Run simulation and generate images in one step
images, systems = simulate_and_image(params, psf;
                                    photons=2000.0,
                                    bg=10.0,
                                    frame_integration=5)
```
"""
function simulate_and_image(params::SmoluchowskiParams, psf::AbstractPSF;
                          photons::Float64=1000.0, bg::Float64=5.0,
                          frame_integration::Int=1, poisson_noise::Bool=true)
    
    # Run diffusion simulation
    systems = simulate(params)
    
    # Generate image sequence
    images = gen_image_sequence(psf, systems, 
                              photons=photons,
                              bg=bg,
                              frame_integration=frame_integration,
                              poisson_noise=poisson_noise)
    
    return images, systems
end