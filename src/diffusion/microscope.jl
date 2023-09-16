# Generating microscope camera data from simulation data 

function gen_image(psf::MicroscopePSFs.PSF,
    states::MoleculeHistory, camera::SMLMSim.Camera, framenum::Int64;
    photons::Float64=1000.0, bg::Float64=5.0, poissonnoise::Bool=true)

    x = [mol.x for mol in states.frames[framenum].molecules]
    y = [mol.y for mol in states.frames[framenum].molecules]

    #convert to pixels
    x = x./camera.pixelsize
    y = y./camera.pixelsize

    points = [(x[i], y[i]) for i in eachindex(x)]
    
    # make empty image
    image = zeros(Float64, camera.ypixels, camera.xpixels)
     
    # calc psf in a 10x10 pixel region around each molecule
    psfbox = 10
    for point in points
        #sub roi
        starty = max(1, Int64(round(point[2]-psfbox)))
        endy = min(camera.ypixels, Int64(round(point[2]+psfbox)))
        startx = max(1, Int64(round(point[1]-psfbox)))
        endx = min(camera.xpixels, Int64(round(point[1]+psfbox)))
        # get the subroi (x,y) pairs. 
        subroi = [(j,i) for i in starty:endy, j in startx:endx]
        # get the psf
        psf_image = photons .* MicroscopePSFs.pdf(psf, subroi, point)
        # add the psf to the image in the correct place
        image[starty:endy,startx:endx] .+= psf_image
    end

    image .+= bg
    if poissonnoise
        for idx in eachindex(image)
            image[idx] = rand(Poisson(image[idx]))
        end
    end
    return image
end

function gen_image_stack(psf::MicroscopePSFs.PSF,
    states::MoleculeHistory, camera::SMLMSim.Camera;
    photons::Float64=1000.0, bg::Float64=5.0, poissonnoise::Bool=true,
    frame_integration::Int64=1)

    # calculate the number of frames in the stack
    num_frames_raw = length(states.frames)
    num_frames = Int64(ceil(num_frames_raw/frame_integration))
    
    # initialize the image stack
    image_stack = zeros(Float64, camera.xpixels, camera.ypixels, num_frames)
    
    # loop over the output frames
    Threads.@threads for i in 1:num_frames
        # loop over the integration frames
        for j in 1:frame_integration
            # calculate the frame number
            frame_num = (i-1)*frame_integration + j
            # check if the frame number is valid
            if frame_num > num_frames_raw
                break
            end
            # calculate the image
            image = gen_image(psf, states, camera, frame_num; photons=photons, bg=bg, poissonnoise=poissonnoise)
            # add the image to the stack
            image_stack[:,:,i] += image
        end
    end

    return image_stack
end

