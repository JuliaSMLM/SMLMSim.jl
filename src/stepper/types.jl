# Parameter types (immutable, user-facing) and runtime state (internal)

# Tolerance on the sum of the mobility fractions
const _FRACTION_TOL = 1e-9

# Transitive closure of the off-diagonal rate graph: is every state reachable from every other?
function _irreducible(q::AbstractMatrix)
    n = size(q, 1)
    reach = [i == j || q[i, j] > 0 for i in 1:n, j in 1:n]
    for k in 1:n, i in 1:n, j in 1:n
        reach[i, k] && reach[k, j] && (reach[i, j] = true)
    end
    return all(reach)
end

# Stationary distribution of an irreducible rate matrix (π q = 0, Σπ = 1)
function _stationary(q::Matrix{Float64})
    n = size(q, 1)
    n == 1 && return [1.0]
    A = Matrix(q')
    A[n, :] .= 1.0
    b = zeros(n)
    b[n] = 1.0
    return A \ b
end

"""
    Population(; name, layer, density, lifetime, birth_rate, mobility, fluor,
               brightness_sigma, budget, multiplicity, z, psf, brightness_jitter, jitter_time)

One kind of emitter in a [`SimWorld`](@ref): in-focus diffusers, immobile clusters,
out-of-focus (OOF) emitters or haze. Built by keyword. Units are μm, s and photons; rates
are per second.

# Keywords
- `name::Symbol = :emitters`: label.
- `layer::Symbol = :signal`: `:signal` or `:oof`, the noise-free map the population renders into.
- `density::Float64`: μm⁻², initial Poisson density over the world box (required).
- `lifetime::Float64 = Inf`: s, mean residence time (exponential departure, independent of excitation).
- `birth_rate::Float64 = isfinite(lifetime) ? density/lifetime : 0.0`: μm⁻² s⁻¹. Holds `density`
  steady only if `budget = Inf`; with a finite budget it settles at `birth_rate` times the mean residence.
- `mobility = [(1.0, 0.0)]`: `(fraction, D μm²/s)` components, one drawn per emitter.
- `fluor::GenericFluor` (required): `γ` is the median photon rate in state 1 at relative intensity 1
  and `q` the CTMC rates at intensity 1. It must be irreducible when it has more than one state.
- `brightness_sigma::Float64 = 0.0`: σ of log γ, per emitter: `γ_i = γ exp(brightness_sigma ξ)`.
- `budget::Float64 = Inf`: mean photons per fluorophore before bleaching (exponential).
- `multiplicity::Int = 1`: fluorophores per emitter (dimming clusters); above 1 needs a one-state `q`.
- `z::Tuple{Float64,Float64} = (0.0, 0.0)`: emitter height ~ U(z[1], z[2]) μm from focus. The
  excitation function reads it.
- `psf` (required): a `GaussianPSF` with fixed σ (μm) or a `StampTable` whose z range covers `z`.
- `brightness_jitter::Float64 = 0.0`, `jitter_time::Float64 = 0.01` (s): frame-to-frame brightness
  fluctuation. Emitter i emits at `γ_i exp(X_i(t))`, where `X_i` is an Ornstein-Uhlenbeck process with
  mean 0, stationary sd `brightness_jitter` and correlation time `jitter_time`. It is drawn stationary at
  birth and held over each sub-step, like motion and excitation. The median rate stays `γ_i`; the mean
  rises by `exp(brightness_jitter^2/2)`. It scales emission only, not the CTMC rates, and the budget
  is spent by emitted photons, so bleaching follows it. With `brightness_jitter = 0` nothing is drawn.
  Real one-molecule movies show a within-track sd of log photons of 0.37-0.48 per 10 ms frame, against
  0.24 without jitter (`PPIDetect/dev/output/t15/sim_vs_real.md`).

# Conventions
`z > 0` points into the sample. A far out-of-focus population uses `layer = :oof`, a blinking
`fluor`, a finite `budget`, and `z = (0.5, 1.0)`. A mirror population at `z = (-1.0, -0.5)` is
meaningful for epi and HILO illumination only. [`UniformExcitation`](@ref) returns 1 at every
z, the baseline for calibrated populations; a z-dependent law is supplied by the caller's
excitation function.

Throws `ArgumentError` for mobility fractions that do not sum to 1, `multiplicity > 1`
with a multi-state `q`, a multi-state `q` that is not irreducible, births with both
`lifetime` and `budget` infinite, a `StampTable` that does not cover `z`, or an unknown `layer`.
"""
struct Population
    name::Symbol
    layer::Symbol
    density::Float64
    lifetime::Float64
    birth_rate::Float64
    mobility::Vector{Tuple{Float64,Float64}}
    fluor::GenericFluor
    brightness_sigma::Float64
    budget::Float64
    multiplicity::Int
    z::Tuple{Float64,Float64}
    psf::Union{GaussianPSF{Float64},StampTable}
    brightness_jitter::Float64
    jitter_time::Float64

    function Population(name, layer, density, lifetime, birth_rate, mobility, fluor,
                        brightness_sigma, budget, multiplicity, z, psf, brightness_jitter, jitter_time)
        layer in (:signal, :oof) || throw(ArgumentError("layer must be :signal or :oof, got :$layer"))
        density >= 0 || throw(ArgumentError("density must be >= 0"))
        lifetime > 0 || throw(ArgumentError("lifetime must be > 0"))
        budget > 0 || throw(ArgumentError("budget must be > 0"))
        birth_rate >= 0 || throw(ArgumentError("birth_rate must be >= 0"))
        brightness_sigma >= 0 || throw(ArgumentError("brightness_sigma must be >= 0"))
        multiplicity >= 1 || throw(ArgumentError("multiplicity must be >= 1"))
        (isfinite(brightness_jitter) && brightness_jitter >= 0) ||
            throw(ArgumentError("brightness_jitter must be finite and >= 0"))
        jitter_time > 0 || throw(ArgumentError("jitter_time must be > 0"))
        z[1] <= z[2] || throw(ArgumentError("z must satisfy z[1] <= z[2]"))
        isempty(mobility) && throw(ArgumentError("mobility must not be empty"))
        for (f, D) in mobility
            (f >= 0 && D >= 0) || throw(ArgumentError("mobility fractions and D must be >= 0"))
        end
        abs(sum(first, mobility) - 1) <= _FRACTION_TOL ||
            throw(ArgumentError("mobility fractions must sum to 1, got $(sum(first, mobility))"))
        q = Matrix{Float64}(fluor.q)
        (size(q, 1) == size(q, 2) && size(q, 1) >= 1) || throw(ArgumentError("q must be a square matrix"))
        for s in 1:size(q, 1)
            off = sum(q[s, j] for j in 1:size(q, 2) if j != s; init=0.0)
            (all(q[s, j] >= 0 for j in 1:size(q, 2) if j != s) && abs(q[s, s] + off) <= 1e-8 * max(1.0, off)) ||
                throw(ArgumentError("row $s of q must have non-negative off-diagonals summing to -q[$s,$s]"))
        end
        nst = size(q, 1)
        multiplicity > 1 && nst > 1 &&
            throw(ArgumentError("multiplicity > 1 needs a one-state q, got $nst states"))
        nst > 1 && !_irreducible(q) &&
            throw(ArgumentError("q must be irreducible (no absorbing or unreachable state); use budget for permanent loss"))
        birth_rate > 0 && !isfinite(lifetime) && !isfinite(budget) &&
            throw(ArgumentError("birth_rate > 0 with infinite lifetime and budget would grow without bound"))
        if psf isa StampTable
            zlo, zhi = first(psf.zs), last(psf.zs)
            (z[1] >= zlo - 1e-12 && z[2] <= zhi + 1e-12) ||
                throw(ArgumentError("the StampTable's z range $zlo..$zhi does not cover z = $z"))
        end
        return new(name, layer, density, lifetime, birth_rate, mobility, fluor,
                   brightness_sigma, budget, multiplicity, z, psf, brightness_jitter, jitter_time)
    end
end

function Population(; name::Symbol=:emitters, layer::Symbol=:signal, density::Real,
                    lifetime::Real=Inf,
                    birth_rate::Real=isfinite(lifetime) ? density / lifetime : 0.0,
                    mobility=[(1.0, 0.0)], fluor::GenericFluor, brightness_sigma::Real=0.0,
                    budget::Real=Inf, multiplicity::Integer=1, z=(0.0, 0.0), psf,
                    brightness_jitter::Real=0.0, jitter_time::Real=0.01)
    mob = Tuple{Float64,Float64}[(Float64(f), Float64(D)) for (f, D) in mobility]
    return Population(name, layer, Float64(density), Float64(lifetime), Float64(birth_rate), mob, fluor,
                      Float64(brightness_sigma), Float64(budget), Int(multiplicity),
                      (Float64(z[1]), Float64(z[2])), psf isa GaussianPSF ? GaussianPSF(Float64(psf.σ)) : psf,
                      Float64(brightness_jitter), Float64(jitter_time))
end

"""
    BackgroundModel(; level, stretch, jitter, feature_size, contrast, correlation_time, illumination_width)

Parameters of the structured background. Built by keyword.

- `level = 0.0`: photons/px/s; a `Real`, or any distribution drawn per stretch with `rand(rng, level)`.
- `stretch::Float64 = Inf`: s; a new level is drawn at `t0 + k stretch`, counted from the world's `t0`, not from the end of `t_burn`.
- `jitter::Float64 = 0.0`: sd of the iid per-exposure multiplier `max(0, 1 + jitter ξ)`.
- `feature_size::Float64 = 0.8`: μm, σ of the pattern's spatial autocorrelation.
- `contrast::Float64 = 0.0`: sd of log pattern; 0 is flat.
- `correlation_time::Float64 = Inf`: s, OU time constant of the pattern; `Inf` is static.
- `illumination_width::Float64 = Inf`: μm, σ of a broad Gaussian profile (mean 1 over the FOV); `Inf` is flat.
"""
Base.@kwdef struct BackgroundModel{L}
    level::L = 0.0
    stretch::Float64 = Inf
    jitter::Float64 = 0.0
    feature_size::Float64 = 0.8
    contrast::Float64 = 0.0
    correlation_time::Float64 = Inf
    illumination_width::Float64 = Inf
end

"""
    UniformExcitation()

Relative excitation intensity 1 at every position, height and time: the baseline for
populations whose `fluor.γ` and `fluor.q` are calibrated at intensity 1.
"""
struct UniformExcitation end

(::UniformExcitation)(x::Float64, y::Float64, z::Float64, t::Float64) = 1.0

"""
    EvanescentExcitation(; depth = 0.1, stray = 0.0)

TIRF excitation: relative intensity `I(z) = stray + (1 - stray) exp(-max(z, 0)/depth)`, 1 at the
glass. `z` is the emitter's height (the [`Population`](@ref) convention, z > 0 into the sample); the
focal plane is at the glass, so z <= 0 gets 1. `depth` (μm, > 0) is the 1/e intensity depth, about
0.08-0.3 μm; `stray` (in [0, 1]) is the fraction of the glass intensity that is propagating
(scattered) light and reaches every height. It applies to every population, `:oof` included, so a
population calibrated under [`UniformExcitation`](@ref) needs its γ rescaled. Like every excitation it
scales emission and the state-1 exit rate.
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

(e::EvanescentExcitation)(x::Float64, y::Float64, z::Float64, t::Float64) =
    e.stray + (1 - e.stray) * exp(-max(z, 0.0) / e.depth)

"""
    next_switch(excitation, t) -> Float64

First time greater than `t` at which `excitation` changes discontinuously, or `Inf`. The
default method returns `Inf`. An excitation is right-continuous: at a switch time it returns
the new value. [`step!`](@ref) evaluates the excitation of each emitter at every declared
switch, so a declared switch acts at its exact time; an undeclared change acts from the
next sub-step start.
"""
next_switch(excitation, t::Float64) = Inf

# Runtime state of the structured background (see BackgroundModel)
mutable struct BackgroundState
    level_src::Any                        # model.level, kept apart: reading it through the parametric model boxes
    stretch::Float64                      # concrete copies of the model's numbers, so reads do not box
    jitter::Float64
    contrast::Float64
    tau::Float64
    t0::Float64
    t_prev::Float64
    stretch_idx::Int
    level::Float64                        # L of the current stretch
    frame_level::Float64                  # L J exposure of the last step
    draw::Vector{Float64}                 # scratch for level draws
    P::Matrix{Float64}                    # illumination, mean 1
    g::Matrix{Float64}                    # lattice node values, N(0, 1)
    tmp::Matrix{Float64}                  # rows mixed, columns still on the lattice
    iy::Vector{Int}
    ix::Vector{Int}
    wy::Matrix{Float64}
    wx::Matrix{Float64}
end

"""
    FrameTruth

One row of [`frame_truth`](@ref): what one emitter did during the last exposure `[t_a, t_b)`,
`T = t_b - t_a`. `frame` is the exposure number, `id` the emitter, `pop` its index into
`world.pops`, `m` the fluorophores left at the end of presence. `x`, `y` (μm) are the
photon-weighted mean position (the presence-weighted mean if no photons, else the reference
position), `z` its height. `photons` is the emitted total, `lit` the fraction of `T` in state 1
with `m > 0`, `excitation` the presence-weighted mean relative intensity (NaN if never present).
`t_birth`, `t_bleach` and `t_depart` are event times inside the exposure, else NaN. The
`partner`, `partner_pop`, `bound`, `t_form`, `t_break`, `lit_bound` and `vis_*` fields are 0,
NaN or false without dimers. `overlap` is set by `SimWorld(...; merge_radius)`.
"""
struct FrameTruth
    frame::Int
    id::Int
    pop::Int32
    m::Int32
    x::Float64
    y::Float64
    z::Float64
    photons::Float64
    lit::Float64
    excitation::Float64
    t_birth::Float64
    t_bleach::Float64
    t_depart::Float64
    partner::Int
    partner_pop::Int32
    bound::Float64
    t_form::Float64
    t_break::Float64
    lit_bound::Float64
    vis_form::Bool
    vis_bound::Bool
    overlap::Bool
end

# Emitters of one Population (struct of arrays, capacity-preallocated, swap-remove).
mutable struct PopState
    p::Population
    γ0::Float64
    q::Matrix{Float64}
    exitrate::Vector{Float64}
    π0::Vector{Float64}
    sigma_px::Float64                     # Gaussian σ in pixels (unused for a stamp table)
    rwin::Int                             # render half-window in pixels (Gaussian)
    stamp::Union{Nothing,StampTable{StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64},Int}}}
    n::Int
    id::Vector{Int}
    x::Vector{Float64}; y::Vector{Float64}; z::Vector{Float64}
    D::Vector{Float64}
    γ::Vector{Float64}
    m::Vector{Int32}
    state::Vector{UInt8}
    clock::Vector{Float64}
    budget::Vector{Float64}
    t_depart::Vector{Float64}
    t_birth::Vector{Float64}
    lj::Vector{Float64}                   # log-brightness multiplier X_i (brightness_jitter)
    xr::Vector{Float64}; yr::Vector{Float64}          # frame accumulators (see R9): reference position,
    sx::Vector{Float64}; sy::Vector{Float64}          # photon-weighted offsets,
    sxp::Vector{Float64}; syp::Vector{Float64}        # presence-weighted offsets,
    sph::Vector{Float64}                              # photons,
    t_present::Vector{Float64}; t_lit::Vector{Float64}; sI::Vector{Float64}
    t_bleach_f::Vector{Float64}
    t_next_birth::Float64
end

"""
    SimWorld(rng, camera, pops; n_sub, background = nothing, boundary = :reflecting, margin, t0 = 0.0)

The runtime state for [`step!`](@ref): the emitters of every [`Population`](@ref), the
camera geometry and the noise-free layers. `rng` drives kinetics only; camera noise is a
separate call on the caller's own RNG.

- `n_sub` (required): sub-steps per exposure. Motion and excitation are held over a
  sub-step, so this sets the blur; event times are exact.
- `background`: `nothing` or a [`BackgroundModel`](@ref), which fills the `structured` layer.
- `boundary`: `:reflecting` or `:periodic`.
- `margin`: μm added around the field of view to form the world box. The default is the
  largest kernel half-width among the populations (5σ for a Gaussian, the stamp radius for a
  stamp table), so light from emitters outside the field of view enters correctly. `margin = 0`
  puts the walls at the field-of-view edge.
- `t0`: start time in s.
- `merge_radius`: μm; when > 0, truth rows of `:signal` emitters closer than this that are not partners
  get `overlap = true`.

Pixels must be uniform and square. Initial emitters are the steady ensemble only for
`multiplicity = 1`, `budget = Inf` and excitation 1; otherwise start the first exposure at
`t0 + t_burn` so the gap advance runs the kinetics unrecorded, with `t_burn` at least 5 mean
residence times.
"""
mutable struct SimWorld{R<:AbstractRNG,C<:AbstractCamera}
    rng::R
    camera::C
    px::Float64
    x0::Float64
    y0::Float64
    box::NTuple{4,Float64}
    boundary::Symbol
    n_sub::Int
    t::Float64
    frame::Int
    pops::Vector{PopState}
    bg::Union{Nothing,BackgroundState}
    signal::Matrix{Float64}
    oof::Matrix{Float64}
    structured::Matrix{Float64}
    expected::Matrix{Float64}
    n_gap::Int
    buf::RenderBuffer
    id_counter::Int
    truth::Vector{FrameTruth}             # rows of the last step!, valid to n_truth
    n_truth::Int
    merge_radius::Float64
    overlap_scratch::Vector{Bool}
    t_a::Float64                          # the recorded exposure
    t_b::Float64
end
