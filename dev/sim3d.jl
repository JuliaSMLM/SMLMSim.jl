using Revise 
using SMLMSim
using SMLMData
using CairoMakie

# Simulation paramters use physical units
# smld structures are in units of pixels and frames 

smld_true, smld_model, smld_noisy = SMLMSim.sim(;
    ρ=1.0,
    σ_PSF=[0.13, 0.13, 0.3], 
    minphotons=50,
    ndatasets=10,
    nframes=1000,
    framerate=50.0, # 1/s
    pattern=SMLMSim.Nmer3D(),
    molecule=SMLMSim.GenericFluor(; q=[0 50; 1e-2 0]), #1/s 
    camera=SMLMSim.IdealCamera(; ypixels=256, xpixels=128, pixelsize=0.1) #pixelsize is microns
)

fig = Figure()
ax = Axis(fig[1, 1],
    xlabel="x",
    ylabel="y",
    aspect=DataAspect())

scatter!(smld_noisy.x, smld_noisy.y; color = smld_noisy.z)
fig
