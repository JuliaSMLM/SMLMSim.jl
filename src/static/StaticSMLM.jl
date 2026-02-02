"""
    StaticSMLM

Module for simulating static (non-diffusing) SMLM data with blinking kinetics.

This module provides functionality for:
1. Generating spatial distributions of emitters based on patterns
2. Simulating fluorophore blinking based on stochastic kinetic models
3. Adding realistic localization uncertainties
4. Creating complete SMLM datasets

# Usage
```julia
using SMLMSim.StaticSMLM
```
"""
module StaticSMLM

using SMLMData
using Distributions

# Import directly from Core module instead of main package
import ..Core: Pattern, Pattern2D, Pattern3D, uniform2D, uniform3D
import ..Core: Molecule, GenericFluor, kinetic_model
import ..Core: Nmer2D, Nmer3D
import ..Core: AbstractSim, SMLMSimParams
import ..Core: AbstractLabeling, FixedLabeling, PoissonLabeling, BinomialLabeling
import ..Core: n_fluorophores, apply_labeling
import ..simulate
import ..SimInfo

include("parameters.jl")
include("coordinate_noise.jl")
include("simulation.jl")

# Export functions and types
export
    # Core types
    StaticSMLMParams,
    
    # Simulation functions
    simulate,
    apply_noise

end # module
