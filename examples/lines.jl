## Simulation of SMLM Data using randomly placed Line2D Pattern

using SMLMSim
using SMLMData
using CairoMakie #remote
# using GLMakie #local


# Simulation parameters use physical units
# smld structures are in units of pixels and frames 

smld_true, smld_model, smld_noisy = SMLMSim.sim(;
    ρ=1.0,
    σ_PSF=0.13, #micron 
    minphotons=50,
    ndatasets=10,
    nframes=1000,
    framerate=50.0, # 1/s
    pattern=SMLMSim.Line2D(),
    molecule=SMLMSim.GenericFluor(; q=[0 50; 1e-2 0]), #1/s 
    camera=SMLMSim.IdealCamera(; xpixels=256, ypixels=256, pixelsize=0.1) #pixelsize is microns
)

fig = Figure()
ax = Axis(fig[1, 1],
    xlabel="x",
    ylabel="y",
    aspect=AxisAspect(1))

scatter!(smld_noisy.x, smld_noisy.y)
fig




