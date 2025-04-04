"""
    Interface Module

This module provides a unified interface for SMLM simulations,
dispatching to the appropriate specialized modules based on the
simulation type and parameters.
"""

"""
    AbstractSim

Abstract type for all simulation types in SMLMSim.
Concrete subtypes should implement their own simulate methods.
"""
abstract type AbstractSim end

"""
    simulate(sim::AbstractSim; kwargs...)

Generic interface for all simulation types.
Dispatches to the appropriate method based on the concrete simulation type.

# Arguments
- `sim::AbstractSim`: The simulation configuration object
- `kwargs...`: Additional keyword arguments specific to the simulation type

# Returns
- The result of the specific simulation method

# Example
```julia
# Create a static SMLM simulation configuration
params = StaticSMLMParams(
    ρ = 1.0,
    σ_psf = 0.13
)
sim_obj = StaticSim(params)

# Run the simulation
results = simulate(sim_obj)
```
"""
function simulate(sim::AbstractSim; kwargs...)
    # This is a generic interface that will dispatch to the
    # appropriate method based on the concrete type of sim
    error("No simulate method implemented for $(typeof(sim))")
end

