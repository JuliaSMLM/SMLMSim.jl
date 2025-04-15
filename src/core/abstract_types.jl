"""
    AbstractSim

Abstract type for all simulation types in SMLMSim.
Concrete subtypes should implement their own simulate methods.
"""
abstract type AbstractSim end

"""
    SMLMSimParams <: AbstractSim

Abstract type for all SMLM simulation parameter types.
Provides a common parent for different types of SMLM simulations.
"""
abstract type SMLMSimParams <: AbstractSim end
