## Simulation of SMLM Data using Nmer Pattern

using Revise
using SMLMSim
using SMLMData
using GLMakie

# Simulation parameters use physical units
# smld structures are in units of pixels and frames 

smld_true, smld_model, smld_noisy=SMLMSim.sim(;
ρ=1.0,
σ_PSF=.13, #micron 
minphotons=50,
ndatasets=10,
nframes=1000,
framerate=50.0, # 1/s
pattern=SMLMSim.Line2D(),
molecule=SMLMSim.GenericFluor(;q=[0 50; 1e-2 0]), #1/s 
camera=SMLMSim.IdealCamera(;xpixels=256,ypixels=256,pixelsize=0.1) #pixelsize is microns
)

plt=Figure()
ax=Axis(plt[1,1])
ax.aspect = AxisAspect(1)
GLMakie.scatter!(smld_noisy.x, smld_noisy.y)
display(plt)

