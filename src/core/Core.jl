"""
    Core

Core module with shared utilities for SMLM simulation.

This module contains fundamental components that can be used across
different simulation types (static, diffusion, etc.) including:

1. Abstract types for simulation
2. Molecule and pattern definitions
3. CTMC (Continuous Time Markov Chain) for stochastic state transitions
4. Photophysics modeling for blinking kinetics and detection

# Usage
```julia
using SMLMSim.Core
```
"""
module Core

using SMLMData
using Distributions
using LinearAlgebra
using ..SMLMSim

include("molecules.jl")
include("patterns.jl")
include("ctmc.jl")
include("photophysics.jl")

# Export abstract types
export
    AbstractSim

# Export molecule types
export
    # Abstract molecule types
    Molecule,
    
    # Concrete molecule types
    GenericFluor

# Export pattern types and functions
export
    # Abstract pattern types
    Pattern,
    Pattern2D,
    Pattern3D,

    # Concrete patterns
    Nmer2D,
    Nmer3D,
    Line2D,
    Line3D,

    # Pattern generation
    uniform2D,
    uniform3D,
    rotate!

# Export CTMC types and functions
export
    CTMC,
    get_state,
    get_next
    
# Export photophysics functions
export
    intensity_trace,
    kinetic_model

end # module
