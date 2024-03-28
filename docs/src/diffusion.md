```@meta
CurrentModule = SMLMSim.InteractionDiffusion
DocTestSetup = quote
    using SMLMSim
end
```

# SMLMSim.InteractionDiffusion
## Overview

This module simulates interacting particles within a box. Monomers are consdiered point particles, and dimers and monomers kept at fixed separation. At each time step of the simulation, the following actions are taken:
- If two free monomers are within the reaction radius, they are linked to form dimers.  
- Pre-existing dimers are broken with a probabilty of $p = k_{\mathrm{off}}dt$
- Monomer position and dimer center-of-mass positions are updated with diffusion simulated as isotropic Brownian motion using their respective diffusion constant.  
- Dimers undergo rotational diffusion. 
- Boundary conditions are applied to the monomer position and dimer center-of-mass

## Tutorial

### Running the Simulation
The main interface to the interaction-diffusion simulator is `smoluchowski()`. All arguments are keyword args and can be changed from default.  

Simply run:

```julia
state_history, args = SMLMSim.smoluchowski()
```

The output is a `MoleculeHistory` structure and `args`, which is a structure that contains all the keyword input arguments that are used in the simulation.  

The `MoleculeHistory` contains the position and monomer/dimer state 
of each particle at each time frame.

The full list of options are in the docstring. 

```@docs
SMLMSim.smoluchowski()
```


### Visualizing the Simulation

To see positions of the particles in a single frame:

```@example 
using SMLMSim
state_history, args = SMLMSim.smoluchowski()
framenum = 100
SMLMSim.show_frame(state_history,framenum,args)
SMLMSim.show_frame(state_history,framenum,args,"diffusion/IDFrame.png") # hide
```

```@raw html
<img src="IDFrame.png" alt="Simulation Frame" width="400"/>
```

Monomers appear as blue dots. Dimers appear as red dots. 


To generate an mp4 movie of all frames:

```@example 
using SMLMSim
state_history, args = SMLMSim.smoluchowski()
SMLMSim.gen_movie(state_history,args; filename="defaultsim.mp4")
```

```@raw html
<video width="320" height="240" controls>
  <source src="defaultsim.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
```


### Simulating Microscope Data 

The `MoleculeHistory` structure returned by `smoluchowski()` can 
be use to simulate data as if it was collected in microscope.  This is implemented by the function `gen_image_stack` and requires a `MicroscopePSFs.PSF` and a `SMLMSim.Camera`.  

A camera has a finite integration time that may be long compared to the dynamics of the diffusion simulation.  The 
The example below shows the simulation of a system and the generation of data stack where each camera image is an integration over `sub_sampling` simulation frames. 


```@example 
using SMLMSim
using MicroscopePSFs
using Images

camera_exposure_time = 0.05 # s
box_size = 10 # microns
sub_sampling = 10
dt = camera_exposure_time/sub_sampling # s

state_history, args = SMLMSim.InteractionDiffusion.smoluchowski(; 
    dt=dt, box_size=box_size);

# Setup a Camera
pixelsize = 0.1
pixels = Int64(round(box_size/pixelsize))
camera = SMLMSim.IdealCamera(; xpixels=pixels, ypixels=pixels, pixelsize=pixelsize)

# Setup a PSF
na = 1.3
wavelength = 0.6 # micron
psf = MicroscopePSFs.Airy2D(na,wavelength,pixelsize)

image_stack = SMLMSim.gen_image_stack(psf, state_history, camera; 
    photons=1000.0, 
    bg=5.0, 
    poissonnoise=true, 
    frame_integration=sub_sampling);

# Look at one frame
im = Gray.(image_stack[:,:,1]./maximum(image_stack[:,:,1]))

save("simimageframe.png",im) # hide
```

```@raw html
<img src="simimageframe.png" alt="Image Frame" width="400"/>
```

### Extracting Dimers

We can extract the dimers into a `MoleculeHistory` structure, which can then be used with the above visualization and image generation tools. 

```@setup dimer
using SMLMSim
using MicroscopePSFs
using Images

camera_exposure_time = 0.05 # s
box_size = 10 # microns
sub_sampling = 10
dt = camera_exposure_time/sub_sampling # s

state_history, args = SMLMSim.InteractionDiffusion.smoluchowski(; 
    dt=dt, box_size=box_size);

# Setup a Camera
pixelsize = 0.1
pixels = Int64(round(box_size/pixelsize))
camera = SMLMSim.IdealCamera(; xpixels=pixels, ypixels=pixels, pixelsize=pixelsize)

# Setup a PSF
na = 1.3
wavelength = 0.6 # micron
psf = MicroscopePSFs.Airy2D(na,wavelength,pixelsize)

image_stack = SMLMSim.gen_image_stack(psf, state_history, camera; 
    photons=1000.0, 
    bg=5.0, 
    poissonnoise=true, 
    frame_integration=sub_sampling);
```

```@example dimer
dimer_history = SMLMSim.get_dimers(state_history)

dimer_stack = SMLMSim.gen_image_stack(psf, dimer_history, camera; 
    photons=1000.0, bg=5.0, poissonnoise=true, frame_integration=10);

framenum = length(state_history.frames)
SMLMSim.show_frame(state_history,framenum,args)
SMLMSim.show_frame(state_history,framenum,args,"DimerFrame.png") # hide

all_image = Gray.(image_stack[:,:,end]./maximum(image_stack[:,:,end]))

dimer_image = Gray.(dimer_stack[:,:,end]./maximum(dimer_stack[:,:,end]))

save("all_image.png",all_image) # hide
save("dimer_image.png",dimer_image) # hide
```

```@raw html
<div style="display: flex; flex-direction: row;">
    <img src="DimerFrame.png" alt="DimerFrame" width="300" style="margin-right: 10px;"/>
    <img src="all_image.png" alt="all_image" width="300" style="margin-right: 10px;"/>
    <img src="dimer_image.png" alt="dimer_image" width="300"/>
</div>
```


## API

```@docs
SMLMSim.InteractionDiffusion
```

```@index
Modules = [InteractionDiffusion]
```

```@autodocs
Modules = [InteractionDiffusion]
```
