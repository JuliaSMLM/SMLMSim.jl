```@meta
CurrentModule = SMLMSim
```

# API Reference

This page provides a comprehensive reference for the types and functions in SMLMSim.

## Core Simulation

```@docs
simulate
sim
```

## Pattern Types

### Abstract Types

```@docs
Pattern
Pattern2D
Pattern3D
```

### 2D Patterns

```@docs
Nmer2D
Line2D
```

### 3D Patterns

```@docs
Nmer3D
Line3D
```

### Pattern Generation

```@docs
uniform2D
uniform3D
rotate!
```

## Molecules and Photophysics

### Types

```@docs
Molecule
GenericFluor
```

### Kinetic Models

```@docs
CTMC
getstate
getnext
intensitytrace
kineticmodel
```

## Localization Uncertainty

```@docs
noise
```

## Interaction-Diffusion

### Types

```@docs
DiffusingMolecule
DiffusingMoleculeSystem
SmoluchowskiParams
```

### Simulation Functions

```@autodocs
Modules = [InteractionDiffusion]
Filter = f -> (f === simulate)
```

### Analysis Functions

```@docs
get_dimers
gen_dimer_images
analyze_dimer_fraction
```

### Microscope Image Generation

```@docs
gen_image
gen_image_sequence
```

### Visualization

```@docs
show_frame
visualize_sequence
visualize_simulation
```

## SMLMData Integration

SMLMSim.jl re-exports several types from SMLMData.jl. Below is a summary of the key types for reference.

### Camera Types

```@docs
AbstractCamera
IdealCamera
```

### Emitter Types

```@docs
AbstractEmitter
Emitter2D
Emitter3D
Emitter2DFit
Emitter3DFit
BasicSMLD
```

## Helper Functions

```@docs
kinetic_model
```

## Index

```@index
```