"""
    SMLMSimParams <: AbstractSMLMConfig

Abstract type for all SMLM simulation parameter types.
Inherits from SMLMData.AbstractSMLMConfig to participate in the
ecosystem-wide (Config, Info, Data) tuple pattern.
"""
abstract type SMLMSimParams <: AbstractSMLMConfig end
