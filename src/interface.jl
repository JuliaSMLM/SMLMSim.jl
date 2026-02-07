"""
    Interface Module

This module provides a unified interface for SMLM simulations,
dispatching to the appropriate specialized modules based on the
simulation type and parameters.
"""

# Import the abstract types directly from Core
import ..Core: SMLMSimParams

# Abstract types are now defined in Core.abstract_types.jl

"""
    simulate(sim::SMLMSimParams; kwargs...)

Generic interface for all simulation types.
Dispatches to the appropriate method based on the concrete simulation type.

# Arguments
- `sim::SMLMSimParams`: The simulation configuration object
- `kwargs...`: Additional keyword arguments specific to the simulation type

# Returns
- The result of the specific simulation method

# Example
```julia
# Create a static SMLM simulation configuration
params = StaticSMLMConfig(
    density = 1.0,        # Changed from ρ to density
    σ_psf = 0.13
)

# Run the simulation
results = simulate(params)
```
"""
function simulate(sim::SMLMSimParams; kwargs...)
    # This is a generic interface that will dispatch to the
    # appropriate method based on the concrete type of sim
    error("No simulate method implemented for $(typeof(sim))")
end
