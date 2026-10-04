```@meta
CurrentModule = SMLMSim
```

# Stepper: closed-loop simulation and background models

The stepper advances a [`SimWorld`](@ref) one camera exposure at a time. Each call to
[`SMLMSim.step!`](@ref) returns the noise-free expected photons per pixel; camera noise is a separate
call on the caller's own random number generator. The same machinery builds background movies
with [`gen_background`](@ref).

All configuration types are built by keyword. Units are μm, s and photons.

## Populations

A [`Population`](@ref) is one kind of emitter: in-focus diffusers, immobile clusters,
out-of-focus (OOF) emitters or haze. `z > 0` points into the sample. A population with
`layer = :oof` goes through the same kinetics but renders into its own `oof` layer.

```julia
using SMLMSim, MicroscopePSFs, Distributions, Random

q = [-50.0 50.0; 1.0 -1.0]
camera = SCMOSCamera(128, 128, 0.1, 1.6; offset = 100.0, gain = 2.0, qe = 0.9)

diffusers = Population(name = :diffusers, density = 0.6,
    mobility = [(0.67, 0.37), (0.33, 0.075)],
    fluor = GenericFluor(; γ = 9500.0, q = q), budget = 5e4,
    psf = GaussianPSF(0.0936))

# far out-of-focus blobs: blinking, bleaching, above the focal plane
blobs = Population(name = :oof, layer = :oof, density = 0.04, lifetime = 0.1,
    mobility = [(1.0, 0.25)], brightness_sigma = 0.3,
    fluor = GenericFluor(; γ = 5000.0, q = q), budget = 2e4,
    z = (0.5, 1.0), psf = GaussianPSF(0.39))
```

A mirror population at `z = (-1.0, -0.5)` is meaningful for epi and HILO illumination only.

## Stepping a world

```julia
world = SimWorld(Random.Xoshiro(1), camera, [diffusers, blobs]; n_sub = 8,
                 background = BackgroundModel(level = 87.7, jitter = 0.025))
noise_rng = Random.Xoshiro(2)
T = 0.01
adu = zeros(128, 128)
for k in 1:100
    expected = SMLMSim.step!(world, (k - 1) * T, k * T)
    scmos_noise!(noise_rng, copyto!(adu, expected), camera)
end
```

`SMLMSim.layers(world)` returns the `signal`, `oof`, `structured` and `expected` maps of the last step
without allocating. The excitation is a callable `(x, y, z, t) -> Float64` given as the fourth
argument of `SMLMSim.step!`; [`UniformExcitation`](@ref) is the default,
[`EvanescentExcitation`](@ref) gives a TIRF falloff with height, and [`next_switch`](@ref)
declares an excitation's discontinuities. A `Population`'s `brightness_jitter` adds frame-to-frame
brightness fluctuation (an Ornstein-Uhlenbeck process in log brightness). `step!` and `layers` are public but not exported, because the
names are generic (SciML and Agents.jl export a `step!`); write `SMLMSim.step!` and `SMLMSim.layers`.

For small `brightness_jitter`, the per-frame sd of an emitter's log photons is about `brightness_jitter √g`,
`g = (n + 2 Σ_{k=1}^{n-1} (n - k) ρ^k)/n²`, `n = n_sub`, `ρ = exp(-T/(n jitter_time))`, `T` the exposure
(`g = 1` at `n_sub = 1`; 0.7385 at `n_sub = 10`, `T = jitter_time = 0.01` s). To reproduce a measured per-frame
sd `j`, set `brightness_jitter = j/√g`. The calibration of `PPIDetect/dev/output/t15/sim_vs_real.md`:
the measured within-track sd of log photons per 10 ms frame is 0.37 (Cell9) and 0.48 (Cell1) against a 0.24
baseline without jitter, so `j = √(target² - 0.24²)` is 0.28 for Cell9 and 0.42 for Cell1, then
`brightness_jitter = j/√g`.

## SpotExcitation

[`SpotExcitation`](@ref) is the excitation of one or more focused spots on a baseline, a callable
`(x, y, z, t)` for `SMLMSim.step!`. Each [`Spot`](@ref) has a centre, a Gaussian `σ` (so the waist is
`2σ`), a `gain`, a Rayleigh range `z_R` and optional `t_on` and `t_off` in seconds on the world clock;
`SpotExcitation` adds the baseline `base`, an evanescent fraction `f_evan` with decay length `d_evan`, and
a lateral `tilt = (tan θ, φ)` that walks the beam axis with depth. Everything is by keyword.

```julia
rig = SpotExcitation(; base = 1.0, spots = [Spot(; x = 6.4, y = 6.4, σ = 0.0934, gain = 50.0, z_R = 0.23)])
world2 = SimWorld(Random.Xoshiro(4), camera, [diffusers]; n_sub = 8)
SMLMSim.step!(world2, 0.0, 0.01, rig)
```

A spot adds `gain (σ/s)² exp(-r²/(2s²))` with `s = σ √(1 + (z/z_R)²)`, so the rig spot at gain 50 falls to
about 4 at `z = 0.75` μm. `z_R = Inf` is a z-independent column. A spot contributes only for
`t_on <= t < t_off`, and [`next_switch`](@ref) returns the next such time, so a spot that switches inside an
exposure acts from its exact time, not from the next sub-step. An emitter of one label (`multiplicity = 1`),
one state, no bleach, `brightness_jitter = 0` and `brightness_sigma = 0` (or read `γ` as that emitter's own `γ_i`), at
the centre of a spot (in focus, or with `z_R = Inf` and `f_evan = 0`) that is on from `t_on` through `t_b`, emits
`γ · gain · (t_b - t_on)` more photons in an exposure the spot turns on in; a blink or a bleach while the spot is on
reduces it.
The baseline `base` is always on. A warning at construction names a spot whose `z_R` is more than 2x from
`π (2σ)² n/λ`.

## Pairs

With `SimWorld(...; dimers = DimerKinetics(...))`, emitters of populations with `binds = true` pair within
and across populations. One [`DimerKinetics`](@ref) governs every pair; a population with `binds = false`
never pairs, and a binding population needs `multiplicity <= 1`. Each sub-step runs in this order:

1. a cell list over the present, unbound, binding emitters;
2. formation: each pair in contact (3D, using each emitter's own `z`, minimum image under `:periodic`)
   forms with probability `1 - exp(-k_on h)`, or on first contact for `k_on = Inf`;
3. splits: every pair whose break time falls before the end of the sub-step dissociates, a pair formed in
   step 2 of the same sub-step included;
4. the usual per-population kinetics, births and motion.

A pair that splits in a sub-step is bound when the sub-step's contacts are tested, so it is first tested
again in the next sub-step: dissociate, diffuse, form. `D_dimer = :min` gives the complex the smaller
member `D`, so a pair with an immobile member does not move and that member is an anchor. Pair placement follows
the diffusion path's rules (`docs/src/diffusion/rules.md`: one distance, orientation from the displacement,
rigid-body reflection, the split margin), in the stepper's rectangular box and in-plane.

**The unbinding rule.** After a split, partners closer than `r_react` (3D) are moved apart in-plane along their axis
to just over `r_react` (an anchor stays); partners already `r_react` or more apart keep their positions. For
partners at equal height, away from the walls, with `d_dimer < r_react`, the chance that they are within `r_react`
again at the next sub-step start is the Gaussian mass of the disk of radius `r_react` around one partner, seen from
distance `r_react`, with per-axis variance `2 (D_a + D_b) h` (`2 D h` for a mover and an anchor): about 6% at
`r_react = 0.03` μm, `D = 0.37` μm²/s and `h = 10` ms, and about half at 0.3 μm, as MicroscopeAdapt measured. With
`d_dimer >= r_react` the disk is seen from `d_dimer`; a height difference `Δz` shrinks its radius to
`sqrt(r_react² - Δz²)`; a wall changes the mass. A contact then forms with probability `1 - exp(-k_on h)`. A finite
`k_on` is a binding rate in place of a capture radius; `k_on = Inf` is the 0.7 contact rule.

**Dark members.** Bleaching changes emission only, never binding. A molecule of a binding population that bleaches
(and an unlabeled one, `multiplicity = 0`) stays present, dark, diffusing and pairing until it departs, so dark
pairs form and `vis_form`/`vis_bound` separate visible pairs. A binding population with births needs a finite
`lifetime`, since its bleached molecules leave only by departure. A `binds = false` molecule is removed when it
bleaches, as without `dimers`. `DimerSim` also keeps bleached molecules.

```julia
dimers = DimerKinetics(k_on = 50.0, r_react = 0.03, k_off = 0.5, D_rot = 1.0, d_dimer = 0.02)
world = SimWorld(Random.Xoshiro(1), camera, [diffusers]; n_sub = 8, dimers = dimers, merge_radius = 0.25)
```

### Recipes

Immobile, non-binding monomers: `binds = false` keeps them out of every pair.

```julia
fixed = Population(name = :fixed, density = 0.2, mobility = [(1.0, 0.0)], binds = false,
    fluor = GenericFluor(; γ = 9500.0, q = zeros(1, 1)), psf = GaussianPSF(0.0936))
```

Immobile binding partners (anchors): `D = 0` with `binds = true`. Movers that bind them pair at the contact
rate, and the complex stays put under `D_dimer = :min`.

```julia
anchors = Population(name = :anchors, density = 0.2, mobility = [(1.0, 0.0)], binds = true,
    fluor = GenericFluor(; γ = 9500.0, q = zeros(1, 1)), psf = GaussianPSF(0.0936))
```

Transient landings: two populations with births and exponential lifetimes, a short-lived one and a
long-lived one. Landings split 82:18 with dwell times 0.19 s and 4.8 s, so the steady densities are in the
ratio `0.82 * 0.19 : 0.18 * 4.8`, about 15:85. The dwell-time data are in
`~/julia_shared_dev/LidkeLab/PPIDetect/dev/output/t15/transient_dwell.md`; a power-law dwell model is a
later addition.

```julia
rate = 0.05                                       # landings per μm² per s, all kinds
fluor = GenericFluor(; γ = 9500.0, q = zeros(1, 1))
short = Population(name = :short, density = 0.0, birth_rate = 0.82 * rate, lifetime = 0.19,
                   fluor = fluor, psf = GaussianPSF(0.0936), mobility = [(1.0, 0.0)])
long  = Population(name = :long, density = 0.0, birth_rate = 0.18 * rate, lifetime = 4.8,
                   fluor = fluor, psf = GaussianPSF(0.0936), mobility = [(1.0, 0.0)])
```

## Truth

[`SMLMSim.frame_truth`](@ref)`(world)` returns one [`FrameTruth`](@ref) row for every emitter present at any
time in the last exposure `[t_a, t_b)`, `T = t_b - t_a`, as a view of the world's buffer that stays valid
until the next `step!`. Gap sub-steps write no rows. Every sub-step is a half-open window `[t0, t1)`: an event
exactly at `t_b` belongs to the next exposure. The fields:

- `frame`, `id`, `pop` (index into `world.pops`) and `m`, the fluorophores left at the end of presence;
- `x`, `y`: the photon-weighted mean position (the presence-weighted mean when the emitter emitted nothing),
  wrapped into the box under `:periodic`; `z` is the emitter's height;
- `photons`: the emitted total; `excitation`: the presence-weighted mean relative intensity, `NaN` if never
  present; `lit` is the fraction of `T` the emitter spent in its emitting state, state 1 with `m > 0`. It is
  defined by state, not by light: an emitter in state 1 that receives no excitation counts as lit and emits
  nothing. Excitation enters its value only through the state-1 exit rate and bleaching. The light received
  shows in `excitation` and `photons`;
- `t_birth`, `t_bleach` and `t_depart`: event times inside the exposure, else `NaN`; `t_bleach` is the time
  `m` reached 0, and a bleach during a gap shows only as `m = 0` in the next row;
- `partner` (an id, 0 when unbound at the end or at removal), `partner_pop`, `bound` (the fraction of `T`
  bound), `t_form` and `t_break` (event times inside the exposure, else `NaN`), `lit_bound` (the fraction of
  `T` bound and emitting);
- `lit`, `bound` and `lit_bound` are fractions of `T` in [0, 1]; where `m`, brightness and `I` are constant over
  the emitter's presence (no bleach in the exposure, no jitter), `photons = m γ_i I lit T`;
- `vis_form`: set in the exposure where the pair formed, true when both members were emitting at the start of
  the forming sub-step; `vis_bound`: true when the emitter was bound at some point in the exposure and it and
  its last partner both have `lit_bound > 0`, so a pair that breaks inside the exposure is still marked;
- `overlap`: with `SimWorld(...; merge_radius)`, set on both rows of every pair of `:signal` rows within
  `merge_radius` μm, unless one is the other's partner at the end of the exposure or its last partner in the
  exposure (a pair that broke in it); two partners that broke and each paired again with another molecule in the
  same exposure can still be marked.

Without `dimers` the pair fields are 0, `NaN` or `false`. [`params_dict`](@ref)`(world)` records the world's
configuration as HDF5-attribute-ready values.

## Background models

The structured background is `S(r) = L J P(r) exp(c g(r) - c^2/2) (t_b - t_a)`:

- a level `L` redrawn every `stretch` seconds (`level` is a number or any distribution),
- a per-exposure jitter `J`,
- a pattern `g` of unit-variance B-spline noise whose spatial autocorrelation has sigma
  `feature_size` μm and whose temporal correlation is `exp(-Δt/correlation_time)`,
- a broad Gaussian illumination profile `P` of mean 1.

Because `exp(c g - c^2/2)` has mean 1 and `P` has mean 1 over the field of view, `L J exposure` is the
expected frame mean of `S` (with `contrast > 0` one frame's spatial mean scatters around it). OOF
populations (model B) add their own layer. [`gen_background`](@ref) returns both maps and the
per-frame level:

```julia
T = 0.01
bg = BackgroundModel(level = LogUniform(2 / T, 100 / T), stretch = 200T, jitter = 0.025,
                     feature_size = 0.8, contrast = 0.3, correlation_time = 0.5,
                     illumination_width = 30.0)
structured, oof, level = gen_background(Random.Xoshiro(3), camera, bg, 1000;
                                        oof = [blobs], frame_time = T, t_burn = 0.5)
bgmovie = structured .+ oof
```

Here `level = LogUniform(2 / T, 100 / T)` with `T = 0.01` means 2 to 100 photons per pixel per
10 ms frame.

## Reference

```@docs
SMLMSim.Stepper
Population
BackgroundModel
SimWorld
DimerKinetics
FrameTruth
frame_truth
params_dict
UniformExcitation
Spot
SpotExcitation
step!
layers
next_switch
gen_background
```
