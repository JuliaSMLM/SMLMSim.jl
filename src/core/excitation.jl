# Relative excitation and brightness fluctuation shared by the simulation paths

"""
    EvanescentExcitation(; depth = 0.1, stray = 0.0)

TIRF excitation: relative intensity `I(z) = stray + (1 - stray) exp(-max(z, 0)/depth)`, 1 at the
glass. `z` is the emitter's height (μm, z > 0 into the sample); z <= 0 gets 1. `depth` (μm, > 0)
is the 1/e intensity depth, about 0.08-0.3 μm; `stray` (in [0, 1]) is the fraction of the glass
intensity that is propagating (scattered) light and reaches every height. Callable as
`e(x, y, z, t)`.

One type serves both simulation paths:
- the diffusion path, `DiffusionSMLMConfig(; excitation)` with heights drawn in `z_range`;
- the stepper, `SMLMSim.step!(world, excitation)`, where `z` is each emitter's height (the
  [`Population`](@ref) convention) and the focal plane is at the glass. There it applies to every
  population, `:oof` included, so a population calibrated under [`UniformExcitation`](@ref) needs
  its γ rescaled; like every excitation it scales emission and the state-1 exit rate.
"""
struct EvanescentExcitation
    depth::Float64
    stray::Float64
    function EvanescentExcitation(depth::Real, stray::Real)
        depth > 0 || throw(ArgumentError("depth must be > 0, got $depth"))
        0 <= stray <= 1 || throw(ArgumentError("stray must be in [0, 1], got $stray"))
        return new(Float64(depth), Float64(stray))
    end
end

EvanescentExcitation(; depth::Real=0.1, stray::Real=0.0) = EvanescentExcitation(depth, stray)

(e::EvanescentExcitation)(x::Real, y::Real, z::Real, t::Real) =
    e.stray + (1 - e.stray) * exp(-max(Float64(z), 0.0) / e.depth)

# Exact Ornstein-Uhlenbeck step of a log-brightness multiplier X (stationary sd s, correlation
# time τ) over h: X ← a·X + b·ξ with ξ ~ N(0, 1). Internal.
_ou_coeffs(s::Float64, τ::Float64, h::Float64) = (exp(-h / τ), s * sqrt(-expm1(-2h / τ)))
