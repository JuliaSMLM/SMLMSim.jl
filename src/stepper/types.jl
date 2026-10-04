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
               brightness_sigma, budget, multiplicity, z, psf, brightness_jitter, jitter_time, binds)

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
  A molecule that never bleaches (`multiplicity = 0` or `fluor.γ = 0`), and in a world with `dimers` any molecule
  of a `binds = true` population (it stays after bleaching), leaves only by departure, so births need a finite
  `lifetime`; the count then settles at `birth_rate × lifetime` per μm² (in a world with `dimers`, a binding
  molecule bound to a longer-lived partner departs with the pair, so its count settles higher; see `DimerKinetics`).
- `mobility = [(1.0, 0.0)]`: `(fraction, D μm²/s)` components, one drawn per emitter.
- `fluor::GenericFluor` (required): `γ` is the median photon rate in state 1 at relative intensity 1
  and `q` the CTMC rates at intensity 1. It must be irreducible when it has more than one state.
- `brightness_sigma::Float64 = 0.0`: σ of log γ, per emitter: `γ_i = γ exp(brightness_sigma ξ)`.
- `budget::Float64 = Inf`: mean photons per fluorophore before bleaching (exponential).
- `multiplicity::Int = 1`: fluorophores per emitter (dimming clusters); above 1 needs a one-state `q`.
  0 is an unlabeled molecule: it never emits and its `t_bleach` is NaN; it still diffuses, departs and, with `dimers`, pairs.
- `z::Tuple{Float64,Float64} = (0.0, 0.0)`: emitter height ~ U(z[1], z[2]) μm from focus. The
  excitation function reads it.
- `psf` (required): a `GaussianPSF` with fixed σ (μm) or a `StampTable` whose z range covers `z`.
- `brightness_jitter::Float64 = 0.0`, `jitter_time::Float64 = 0.01` (s): frame-to-frame brightness
  fluctuation. Emitter i emits at `γ_i exp(X_i(t))`, where `X_i` is an Ornstein-Uhlenbeck process with
  mean 0, stationary sd `brightness_jitter` and correlation time `jitter_time`. It is drawn stationary at
  birth and held over each sub-step, like motion and excitation. The median rate stays `γ_i`; the mean
  rises by `exp(brightness_jitter^2/2)`. It scales emission only, not the CTMC rates, and the budget
  is spent by emitted photons, so bleaching follows it. With `brightness_jitter = 0` nothing is drawn.
  For small `brightness_jitter`, the per-frame sd of an emitter's log photons is about `brightness_jitter·√g`,
  `g = (n + 2 Σ_{k=1}^{n−1} (n − k) ρ^k)/n²`, `n = n_sub`, `ρ = exp(−T/(n·jitter_time))`, `T` the exposure (g = 1 at
  `n_sub = 1`; 0.7385 at `n_sub = 10`, `T = jitter_time = 0.01` s). To reproduce a measured per-frame sd `j`, set
  `brightness_jitter = j/√g`.
  Real one-molecule movies show a within-track sd of log photons of 0.37-0.48 per 10 ms frame, against
  0.24 without jitter (`PPIDetect/dev/output/t15/sim_vs_real.md`).
- `binds::Bool = true`: with `SimWorld(...; dimers)`, pairs with other binding emitters, within and across
  populations; `false` never pairs. A binding population needs `multiplicity <= 1` when the world has `dimers`.

# Conventions
`z > 0` points into the sample. A far out-of-focus population uses `layer = :oof`, a blinking
`fluor`, a finite `budget`, and `z = (0.5, 1.0)`. A mirror population at `z = (-1.0, -0.5)` is
meaningful for epi and HILO illumination only. [`UniformExcitation`](@ref) returns 1 at every
z, the baseline for calibrated populations; a z-dependent law is supplied by the caller's
excitation function.

Throws `ArgumentError` for mobility fractions that do not sum to 1, `multiplicity > 1`
with a multi-state `q`, a multi-state `q` that is not irreducible or has more than 255 states, births with an
infinite `lifetime` and a molecule that cannot bleach (`budget` infinite, `multiplicity = 0` or `fluor.γ = 0`), a `StampTable` that does not cover `z`, an unknown `layer`, or a nonfinite
`density`, `birth_rate`, `brightness_sigma`, mobility `D`, `z` or `fluor.γ` (a nonfinite `q` entry fails its row check;
`lifetime`, `budget` and `jitter_time` may be `Inf`).
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
    binds::Bool

    function Population(name, layer, density, lifetime, birth_rate, mobility, fluor,
                        brightness_sigma, budget, multiplicity, z, psf, brightness_jitter, jitter_time, binds)
        layer in (:signal, :oof) || throw(ArgumentError("layer must be :signal or :oof, got :$layer"))
        (density >= 0 && isfinite(density)) || throw(ArgumentError("density must be finite and >= 0"))
        lifetime > 0 || throw(ArgumentError("lifetime must be > 0"))
        budget > 0 || throw(ArgumentError("budget must be > 0"))
        (birth_rate >= 0 && isfinite(birth_rate)) || throw(ArgumentError("birth_rate must be finite and >= 0"))
        (brightness_sigma >= 0 && isfinite(brightness_sigma)) ||
            throw(ArgumentError("brightness_sigma must be finite and >= 0"))
        multiplicity >= 0 || throw(ArgumentError("multiplicity must be >= 0"))
        (isfinite(brightness_jitter) && brightness_jitter >= 0) ||
            throw(ArgumentError("brightness_jitter must be finite and >= 0"))
        jitter_time > 0 || throw(ArgumentError("jitter_time must be > 0"))
        (isfinite(z[1]) && isfinite(z[2]) && z[1] <= z[2]) ||
            throw(ArgumentError("z must be finite and satisfy z[1] <= z[2]"))
        (isfinite(fluor.γ) && fluor.γ >= 0) || throw(ArgumentError("fluor.γ must be finite and >= 0"))
        isempty(mobility) && throw(ArgumentError("mobility must not be empty"))
        for (f, D) in mobility
            (f >= 0 && D >= 0 && isfinite(D)) ||
                throw(ArgumentError("mobility fractions must be >= 0 and D finite and >= 0"))
        end
        abs(sum(first, mobility) - 1) <= _FRACTION_TOL ||
            throw(ArgumentError("mobility fractions must sum to 1, got $(sum(first, mobility))"))
        q = Matrix{Float64}(fluor.q)
        all(isfinite, q) || throw(ArgumentError("q must be finite"))
        (size(q, 1) == size(q, 2) && size(q, 1) >= 1) || throw(ArgumentError("q must be a square matrix"))
        for s in 1:size(q, 1)
            off = sum(q[s, j] for j in 1:size(q, 2) if j != s; init=0.0)
            (all(q[s, j] >= 0 for j in 1:size(q, 2) if j != s) && abs(q[s, s] + off) <= 1e-8 * max(1.0, off)) ||
                throw(ArgumentError("row $s of q must have non-negative off-diagonals summing to -q[$s,$s]"))
        end
        nst = size(q, 1)
        nst <= 255 || throw(ArgumentError("q has $nst states; the CTMC state is stored in a UInt8, so at most 255"))
        multiplicity > 1 && nst > 1 &&
            throw(ArgumentError("multiplicity > 1 needs a one-state q, got $nst states"))
        nst > 1 && !_irreducible(q) &&
            throw(ArgumentError("q must be irreducible (no absorbing or unreachable state); use budget for permanent loss"))
        birth_rate > 0 && !isfinite(lifetime) && !(isfinite(budget) && multiplicity > 0 && fluor.γ > 0) &&
            throw(ArgumentError("birth_rate > 0 with an infinite lifetime needs a molecule that can bleach (finite budget, multiplicity > 0 and fluor.γ > 0); otherwise the count would grow without bound"))
        if psf isa StampTable
            zlo, zhi = first(psf.zs), last(psf.zs)
            (z[1] >= zlo - 1e-12 && z[2] <= zhi + 1e-12) ||
                throw(ArgumentError("the StampTable's z range $zlo..$zhi does not cover z = $z"))
        end
        return new(name, layer, density, lifetime, birth_rate, mobility, fluor,
                   brightness_sigma, budget, multiplicity, z, psf, brightness_jitter, jitter_time, binds)
    end
end

function Population(; name::Symbol=:emitters, layer::Symbol=:signal, density::Real,
                    lifetime::Real=Inf,
                    birth_rate::Real=isfinite(lifetime) ? density / lifetime : 0.0,
                    mobility=[(1.0, 0.0)], fluor::GenericFluor, brightness_sigma::Real=0.0,
                    budget::Real=Inf, multiplicity::Integer=1, z=(0.0, 0.0), psf,
                    brightness_jitter::Real=0.0, jitter_time::Real=0.01, binds::Bool=true)
    mob = Tuple{Float64,Float64}[(Float64(f), Float64(D)) for (f, D) in mobility]
    return Population(name, layer, Float64(density), Float64(lifetime), Float64(birth_rate), mob, fluor,
                      Float64(brightness_sigma), Float64(budget), Int(multiplicity),
                      (Float64(z[1]), Float64(z[2])), psf isa GaussianPSF ? GaussianPSF(Float64(psf.σ)) : psf,
                      Float64(brightness_jitter), Float64(jitter_time), binds)
end

"""
    DimerKinetics(; k_on, r_react, k_off, D_rot, d_dimer, D_dimer = :min)
    DimerKinetics(cfg::DiffusionSMLMConfig; k_on, D_dimer = cfg.diff_dimer)

The pairing kinetics of a [`SimWorld`](@ref) built with `dimers = DimerKinetics(...)`. One instance governs
every pair, within and across populations; a [`Population`](@ref) with `binds = false` never pairs. Built by
keyword. Units are μm, s and rad.

- `k_on`: s⁻¹, the Doi rate at which a pair within `r_react` reacts; `Inf` reacts on first contact
  (the 0.7 contact rule). Must be > 0.
- `r_react`: μm, the 3D contact distance sqrt(Δx² + Δy² + Δz²) using each emitter's own z. Contacts are tested
  at sub-step starts.
- `k_off`: s⁻¹, dissociation rate (0 = never). Measured EGFR dimers dissociate at 0.12-0.27/s (EGF-bound) and
  0.31-1.24/s (unliganded) (Low-Nam 2011, Valley 2015); 10/s is the fast regime simulated by Pryor 2013, not a
  measurement. The dimer diffuses about 6 times slower than the monomer.
- `D_dimer`: μm²/s of the complex centre; `:min` is min(D_i, D_j), so a pair with an immobile member does not move,
  and that member keeps its position (the anchor), or a fixed value.
- `D_rot`: rad²/s, rotational diffusion of the pair axis.
- `d_dimer`: μm, the in-plane (x, y) separation of the members while bound; each keeps its own height, so their 3D
  distance is `sqrt(d_dimer² + Δz²)`.

A pair leaves as one: its members share one departure time drawn from the longer of their two lifetimes (`Inf`
if either is `Inf`), so the shorter-lived member of a cross pair cannot depart while bound and its population's
steady count exceeds `birth_rate * lifetime` by its time bound to longer-lived partners. After a split each
member draws a fresh departure time. Bleaching changes emission only, never binding: a molecule of a `binds = true`
population that bleaches (and an unlabeled one, `multiplicity = 0`) stays present, dark, diffusing and pairing until
it departs. A `binds = false` molecule is removed when it bleaches, as without `dimers`.

After a split, partners closer than `r_react` (3D) are moved apart in-plane along their axis to just over `r_react`
(an anchor stays); partners already `r_react` or more apart keep their positions. For partners at equal height, away
from the walls, with `d_dimer < r_react`, the chance that they are within `r_react` again at the next sub-step start
is the Gaussian mass of the disk of radius `r_react` around one partner, seen from distance `r_react`, with per-axis
variance σ² = 2(D_a + D_b)h (2Dh for a mover and an anchor): about 6% at `r_react` = 0.03 μm, D = 0.37 μm²/s,
h = 10 ms, and about half at 0.3 μm, as MicroscopeAdapt measured. With `d_dimer >= r_react` the disk is seen from
`d_dimer`; a height difference Δz shrinks its radius to `sqrt(r_react² - Δz²)`; a wall changes the mass.
Each contact then forms with probability 1 - exp(-k_on h). A finite `k_on` is
MicroscopeAdapt's "binding rate in place of a capture radius" option; `k_on = Inf` is the 0.7 contact rule.

Throws `ArgumentError` unless `k_on > 0`, `0 < r_react < Inf`, `0 <= k_off < Inf`, `0 <= D_rot < Inf`,
`0 <= d_dimer < Inf` and `D_dimer` is `:min` or a finite real >= 0.
"""
struct DimerKinetics
    k_on::Float64
    r_react::Float64
    k_off::Float64
    D_dimer::Union{Symbol,Float64}
    D_rot::Float64
    d_dimer::Float64

    function DimerKinetics(k_on, r_react, k_off, D_dimer, D_rot, d_dimer)
        k_on, r_react, k_off, D_rot, d_dimer = Float64(k_on), Float64(r_react), Float64(k_off), Float64(D_rot), Float64(d_dimer)
        (D_dimer === :min || D_dimer isa Real) ||
            throw(ArgumentError("D_dimer must be :min or a finite real >= 0, got $D_dimer"))
        D_dimer === :min || (D_dimer = Float64(D_dimer))
        k_on > 0 || throw(ArgumentError("k_on must be > 0 (Inf allowed), got $k_on"))
        (r_react > 0 && isfinite(r_react)) || throw(ArgumentError("r_react must be in (0, Inf), got $r_react"))
        (k_off >= 0 && isfinite(k_off)) || throw(ArgumentError("k_off must be in [0, Inf), got $k_off"))
        (D_rot >= 0 && isfinite(D_rot)) || throw(ArgumentError("D_rot must be in [0, Inf), got $D_rot"))
        (d_dimer >= 0 && isfinite(d_dimer)) || throw(ArgumentError("d_dimer must be in [0, Inf), got $d_dimer"))
        (D_dimer === :min || (D_dimer >= 0 && isfinite(D_dimer))) ||
            throw(ArgumentError("D_dimer must be :min or a finite real >= 0, got $D_dimer"))
        return new(k_on, r_react, k_off, D_dimer, D_rot, d_dimer)
    end
end

DimerKinetics(; k_on::Real, r_react::Real, k_off::Real, D_rot::Real, d_dimer::Real, D_dimer=:min) =
    DimerKinetics(k_on, r_react, k_off, D_dimer, D_rot, d_dimer)

DimerKinetics(cfg::DiffusionSMLMConfig; k_on::Real, D_dimer=cfg.diff_dimer) =
    DimerKinetics(; k_on, r_react=cfg.r_react, k_off=cfg.k_off, D_dimer, D_rot=cfg.diff_dimer_rot,
                  d_dimer=cfg.d_dimer)

"""
    BackgroundModel(; level, stretch, jitter, feature_size, contrast, correlation_time, illumination_width)

Parameters of the structured background. Built by keyword.

- `level = 0.0`: photons/px/s; a `Real`, or any distribution drawn per stretch with `rand(rng, level)`.
- `stretch::Float64 = Inf`: s; a new level is drawn at `t0 + k stretch`, counted from the world's `t0`, not from the end of `t_burn`; every boundary crossed draws a level, including those in a gap or `t_burn`.
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

# EvanescentExcitation is defined in Core (src/core/excitation.jl) and shared with the diffusion path

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
    model::Any                            # the BackgroundModel, kept for params_dict
end

"""
    FrameTruth

One row of [`frame_truth`](@ref): what one emitter did during the last exposure `[t_a, t_b)`,
`T = t_b - t_a`. `frame` is the exposure number, `id` the emitter, `pop` its index into
`world.pops`, `m` the fluorophores left at the end of presence. `x`, `y` (μm) are the
photon-weighted mean position (the presence-weighted mean if no photons, else the reference
position), `z` its height. `photons` is the emitted total, `excitation` the presence-weighted mean relative
intensity (NaN if never present). `lit` is the fraction of `T` the emitter spent in its emitting state, state 1 with
`m > 0`. It is defined by state, not by light: an emitter in state 1 that receives no excitation counts as lit and
emits nothing. Excitation enters its value only through the state-1 exit rate and bleaching. The light received
shows in `excitation` and `photons`. `lit`, `bound` and `lit_bound` are fractions of `T` in [0, 1]. For an emitter whose `m`, brightness and excitation `I`
are constant over its presence (no bleach in the exposure, `brightness_jitter = 0`), `photons = m γ_i I lit T`; a
bleach inside the exposure leaves `photons > 0` with `m = 0`.
`t_birth`, `t_bleach` and `t_depart` are event times inside the exposure, else NaN. The
`partner`, `partner_pop`, `bound`, `t_form`, `t_break`, `lit_bound` and `vis_*` fields are 0,
NaN or false without dimers. `partner` is the partner at the end of the exposure (or at removal),
0 if unbound then. `bound` is the fraction of `T` the emitter is bound, and `lit_bound`
the fraction of `T` it is lit (state 1 with `m > 0`) while bound, so a member lit for 5 ms of a
10 ms exposure while bound has `lit_bound = 0.5`. `vis_form` is true in the exposure where the
pair formed (`t_form` not NaN) when both members were lit at the start of the forming sub-step.
`vis_bound` is true when the emitter was bound at some point in this exposure and it and that
partner (the last one, if several) both have `lit_bound > 0`, so a pair that breaks within the
exposure is marked in it while its rows' `partner` is already 0. `overlap` is set by `SimWorld(...; merge_radius)`.
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
    partner::Vector{Int32}                            # dimers: index within the partner's population, 0 = unbound
    partner_pop::Vector{Int32}                        # index into world.pops of the partner's population
    partner_id::Vector{Int}                           # id of the partner (truth; outlives the link at removal)
    t_form::Vector{Float64}; t_break_due::Vector{Float64}   # -Inf if never bound
    t_bound::Vector{Float64}; t_break_f::Vector{Float64}    # frame accumulators
    t_litb::Vector{Float64}                           # frame accumulator: time bound and emitting
    vis_form::Vector{Bool}                            # both members emitting at the start of the forming sub-step
    partner_f::Vector{Int}                            # id of the last partner (set at formation, 0 from birth; vis pass)
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
- `dimers`: `nothing` or a [`DimerKinetics`](@ref): emitters of `binds = true` populations pair within and across
  populations. Every binding population needs `multiplicity <= 1` and, when it has births, a finite `lifetime` (its
  bleached molecules stay until they depart), and every box side must exceed twice the larger of `d_dimer` and the
  split separation (`r_react` plus a rounding margin), else `ArgumentError`.
- `merge_radius`: μm; when > 0, truth rows of `:signal` emitters closer than this get `overlap = true`,
  unless one is the other's partner at the end of this exposure, or its last partner in this exposure (as for
  a pair that broke in it; two partners that broke and each paired again with another in the same exposure can
  still be marked).

Pixels must be uniform and square. Initial emitters are the steady ensemble only when
births balance departures (`birth_rate * lifetime == density`, as the default `birth_rate` gives, or
`lifetime = Inf` with no births), `multiplicity <= 1`, `budget = Inf`, the excitation is 1 and the world has no
`dimers` (pairs start unformed); otherwise start the first exposure at
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
    vis_partner::Vector{Int}              # per truth row: the emitter's last partner, read when bound in the exposure
    t_a::Float64                          # the recorded exposure
    t_b::Float64
    dimers::Union{Nothing,DimerKinetics}
    head::Vector{Int32}                   # cell list over unbound binding emitters: first entry of each cell
    next::Vector{Int32}                   # per entry: next entry of the same cell
    ent_pop::Vector{Int32}; ent_idx::Vector{Int32}
    ncx::Int; ncy::Int                    # cells per axis
end
