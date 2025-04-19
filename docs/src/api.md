```@meta
CurrentModule = SMLMSim
```

# API Reference

This page provides a comprehensive reference for the types and functions in SMLMSim.

## Core Simulation

```@docs
simulate
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
get_state
get_next
intensity_trace
kinetic_model
```

## Localization Uncertainty

```@docs
noise
apply_noise
```

## Static Simulation

```@docs
StaticSMLMParams
```

## Interaction-Diffusion

### Types

```@docs
DiffusionSMLMParams
DiffusingEmitter2D
DiffusingEmitter3D
```

### Analysis Functions

```@docs
get_dimers
get_monomers
analyze_dimer_fraction
analyze_dimer_lifetime
```

### Microscope Image Generation

```@docs
gen_images
gen_image
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