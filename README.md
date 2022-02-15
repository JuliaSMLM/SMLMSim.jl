# SMLMSim

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/SMLMSim.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/SMLMSim.jl/dev)
[![Build Status](https://github.com/JuliaSMLM/SMLMSim.jl/workflows/CI/badge.svg)](https://github.com/JuliaSMLM/SMLMSim.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaSMLM/SMLMSim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/SMLMSim.jl)

Generate simulated SMLM coordinate data.  Patterns, Cameras, Fluorphores, and SMLD data organization can be configured. 

Simulation parameters use physical units. Resulting smld structures are in units of pixels and frames. 


The high level interface is the `SMLMSim.sim()` function.   

```
using SMLMSim
using SMLMData
using PlotlyJS

smld_true, smld_model, smld_noisy=SMLMSim.sim(;
ρ=1.0, # pattern density (1/micron^2)
σ_PSF=.13, # (micron) 
minphotons=50,
ndatasets=10,
nframes=1000,
framerate=50.0, # (1/s)
pattern=SMLMSim.Nmer2D(),
molecule=SMLMSim.GenericFluor(;q=[0 50; 1e-2 0]), # (1/s) 
camera=SMLMSim.IdealCamera(;xpixels=256,ypixels=256,pixelsize=0.1) # pixelsize is microns
)

plt=PlotlyJS.plot(scattergl(x=smld_noisy.x, y=smld_noisy.y, mode="markers"))
display(plt)
```
