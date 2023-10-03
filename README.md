# SMLMSim

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSMLM.github.io/SMLMSim.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSMLM.github.io/SMLMSim.jl/dev)
[![Build Status](https://github.com/JuliaSMLM/SMLMSim.jl/workflows/CI/badge.svg)](https://github.com/JuliaSMLM/SMLMSim.jl/actions)
[![Coverage](https://codecov.io/gh/JuliaSMLM/SMLMSim.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaSMLM/SMLMSim.jl)

Generate simulated SMLM coordinate data.  Patterns, Cameras, Fluorophores, and SMLD data organization can be configured. 

Simulation parameters use physical units. Resulting `SMLMData.SMLD` structures are in units of pixels and frames. 


The high level interface for simulating SMLM super-resolution coordinate data is the `SMLMSim.sim()` function.   

```
using SMLMSim
using SMLMData

smld_true, smld_model, smld_noisy = SMLMSim.sim(;
    ρ=1.0,
    σ_PSF=0.13, 
    minphotons=50,
    ndatasets=10,
    nframes=1000,
    framerate=50.0, # 1/s
    pattern=SMLMSim.Nmer2D(),
    molecule=SMLMSim.GenericFluor(; q=[0 50; 1e-2 0]), #1/s 
    camera=SMLMSim.IdealCamera(; ypixels=256, xpixels=128, pixelsize=0.1) #pixelsize is microns
)
```
