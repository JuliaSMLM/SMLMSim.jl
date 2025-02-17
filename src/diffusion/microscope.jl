
"""
    gen_image(psf::PSF, system::DiffusingMoleculeSystem, frame::Int;
              poisson_noise::Bool=true)

Generate a microscope image from a diffusing molecule system.

# Arguments
- `psf::PSF`: Point spread function model (e.g., Airy2D, Gaussian2D)
- `system::DiffusingMoleculeSystem`: The molecular system state
- `frame::Int`: Frame number to generate

# Optional Arguments
- `poisson_noise::Bool=true`: Whether to add Poisson noise to final image

# Returns
- Array{Float64,2}: Generated image with dimensions matching system camera
"""
function gen_image(psf::PSF, system::DiffusingMoleculeSystem, frame::Int;
                  poisson_noise::Bool=true)
    
    # Get camera dimensions from pixel edges
    ny = length(system.camera.pixel_edges_y) - 1
    nx = length(system.camera.pixel_edges_x) - 1
    image = zeros(Float64, ny, nx)
    
    # Process each molecule
    for mol in system.molecules
        # Skip if molecule isn't emitting (e.g. part of a non-fluorescent complex)
        mol.photons <= 0 && continue
        
        # Integrate PSF over pixels for this emitter
        pixel_intensities = integrate_pixels(psf, system.camera, mol.emitter)
        
        # Add to image
        image .+= pixel_intensities
    end
    
    # Get background from metadata or use default
    bg = get(system.metadata, "background", 5.0)
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
    gen_image_sequence(psf::PSF, systems::Vector{DiffusingMoleculeSystem};
                      frame_integration::Int=1, poisson_noise::Bool=true)

Generate a sequence of microscope images from diffusing molecule system states.

# Arguments
- `psf::PSF`: Point spread function model
- `systems::Vector{DiffusingMoleculeSystem}`: Sequence of system states

# Optional Arguments
- `frame_integration::Int=1`: Number of simulation frames to integrate per output frame
- `poisson_noise::Bool=true`: Whether to add Poisson noise

# Returns
- Array{Float64,3}: Stack of generated images [ny, nx, frames]
"""
function gen_image_sequence(psf::PSF, systems::Vector{DiffusingMoleculeSystem};
                          frame_integration::Int=1, poisson_noise::Bool=true)
    
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
                mol.emitter.photons *= photon_scale
            end
            
            # Generate and add frame without noise
            frame .+= gen_image(psf, system, idx, poisson_noise=false)
            
            # Restore original photon counts
            for mol in system.molecules
                mol.emitter.photons /= photon_scale
            end
        end
        
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
    simulate_and_image(params::SmoluchowskiParams, psf::PSF;
                      frame_integration::Int=1, poisson_noise::Bool=true)

Run a complete diffusion simulation and generate microscope images.

# Arguments
- `params::SmoluchowskiParams`: Simulation parameters
- `psf::PSF`: Point spread function model

# Optional Arguments
- `frame_integration::Int=1`: Number of simulation frames to integrate per output frame
- `poisson_noise::Bool=true`: Whether to add Poisson noise

# Returns
- Array{Float64,3}: Stack of generated images [ny, nx, frames]
- Vector{DiffusingMoleculeSystem}: System states at each timestep
"""
function simulate_and_image(params::SmoluchowskiParams, psf::PSF;
                          frame_integration::Int=1, poisson_noise::Bool=true)
    
    # Run diffusion simulation
    systems = simulate(params)
    
    # Generate image sequence
    images = gen_image_sequence(psf, systems, 
                              frame_integration=frame_integration,
                              poisson_noise=poisson_noise)
    
    return images, systems
end