# High Level Interface

function sim(;
    ρ=1.0,
    σ_PSF=0.13,
    minphotons=50,
    ndatasets=10,
    nframes=1000,
    framerate=50.0,
    pattern::Pattern=Nmer2D(),
    molecule::Molecule=GenericFluor(; q=[0 50; 1e-2 0]),
    camera::Camera=IdealCamera()
)

    # Simulation sequence
    smld_true = SMLMSim.uniform2D(ρ, pattern, camera.xpixels * camera.pixelsize, camera.ypixels * camera.pixelsize)
    smld_model = SMLMSim.kineticmodel(smld_true, molecule, nframes, framerate; ndatasets, minphotons)
    smld_noisy = SMLMSim.noise(smld_model, σ_PSF)

    return smld_true, smld_model, smld_noisy
end


