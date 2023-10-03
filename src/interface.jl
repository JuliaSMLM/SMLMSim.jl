# High Level Interface

"""
    sim(;
        ρ=1.0,
        σ_PSF=.13,
        minphotons=50,
        ndatasets=10,
        nframes=1000,
        framerate=50.0, 
        pattern::Pattern=Nmer2D(),
        molecule::Molecule=GenericFluor(;q=[0 50; 1e-2 0]),
        camera::Camera=IdealCamera(),
        zrange::Vector{<:Real}=[-1.0, 1.0]
    )

Generate SMLD using simulation parmeters.      
"""
function sim(;
    ρ=1.0,
    σ_PSF=0.13,
    minphotons=50,
    ndatasets=10,
    nframes=1000,
    framerate=50.0,
    pattern::Pattern=Nmer2D(),
    molecule::Molecule=GenericFluor(; q=[0 50; 1e-2 0]),
    camera::Camera=IdealCamera(),
    zrange::Vector{<:Real}=[-1.0, 1.0]
)

    # Simulation sequence

    # check if pattern is 2D or 3D 

    if pattern isa Pattern2D
        coords = SMLMSim.uniform2D(ρ, pattern, camera.xpixels * camera.pixelsize, camera.ypixels * camera.pixelsize)
        smld_true = SMLMData.SMLD2D(camera, coords[1], coords[2], ones(length(coords[1])),
            ones(Int, length(coords[1])), ones(Int, length(coords[1])), Vector(1:length(coords[1])))
    elseif pattern isa Pattern3D
        coords = SMLMSim.uniform3D(ρ, pattern, camera.xpixels * camera.pixelsize, camera.ypixels * camera.pixelsize,
            zrange=zrange)
        smld_true = SMLMData.SMLD3D(camera, coords[1], coords[2], coords[3], ones(length(coords[1])),
            ones(Int, length(coords[1])), ones(Int, length(coords[1])), Vector(1:length(coords[1])))
    else
        error("pattern must be 2D or 3D")
    end


    out = SMLMSim.kineticmodel(coords..., molecule, nframes, framerate; ndatasets, minphotons)
    
    if pattern isa Pattern2D
        smld_model = SMLMData.SMLD2D(camera, out...)
        σ_PSF = σ_PSF / camera.pixelsize
    elseif pattern isa Pattern3D
        smld_model = SMLMData.SMLD3D(camera, out...)
        σ_PSF = σ_PSF .* [1.0 / camera.pixelsize, 1.0 / camera.pixelsize, 1.0]

    end
    
    smld_noisy = SMLMSim.noise(smld_model, σ_PSF)

    return smld_true, smld_model, smld_noisy
end

