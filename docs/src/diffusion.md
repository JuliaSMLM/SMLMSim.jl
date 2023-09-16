```@meta
CurrentModule = SMLMSim.InteractionDiffusion
DocTestSetup = quote
    using SMLMSim
end
```

# SMLMSim.InteractionDiffusion
## Overview

This module simulates interacting particles within a box. At each time step of the simulation, the following actions are taken:
- Chemical species are updated.  If two free monomers are within the reaction radius, they are linked to form dimers.  Each pre-existing dimer has a probabilty of being broken into monomers given by $p = k_{\mathrm{off}}dt$
- Positions are updated. Monomer position and dimer center-of-mass position are updated with diffusion simulated as isotropic Brownian motion using thier respective diffusion constant.  Dimers additionly undergo rotational diffusion. 
- Boundary conditions are applied to the Monomer position and dimer center-of-mass

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
SMLMSim.show_frame(state_history,framenum,args,"IDFrame.png") # hide
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




## API

```@index
```

```@autodocs
Modules = [SMLMSim.InteractionDiffusion]
Order   = [:type, :function]
```