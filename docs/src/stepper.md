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
argument of `SMLMSim.step!`; [`UniformExcitation`](@ref) is the default and [`next_switch`](@ref)
declares its discontinuities. `step!` and `layers` are public but not exported, because the
names are generic (SciML and Agents.jl export a `step!`); write `SMLMSim.step!` and `SMLMSim.layers`.

## Background models

The structured background is `S(r) = L J P(r) exp(c g(r) - c^2/2) (t_b - t_a)`:

- a level `L` redrawn every `stretch` seconds (`level` is a number or any distribution),
- a per-exposure jitter `J`,
- a pattern `g` of unit-variance B-spline noise whose spatial autocorrelation has sigma
  `feature_size` μm and whose temporal correlation is `exp(-Δt/correlation_time)`,
- a broad Gaussian illumination profile `P` of mean 1.

Because `exp(c g - c^2/2)` has mean 1, the frame mean of `S` is `L J exposure`. OOF
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
UniformExcitation
step!
layers
next_switch
gen_background
```
