```@meta
CurrentModule = SMLMSim
```

# Diffusion-Interaction Simulation

## Overview

The Diffusion-Interaction module simulates dynamic molecular processes including diffusion and interactions between particles in a controlled environment. This allows you to model realistic molecular behaviors such as:

- Free diffusion of monomers
- Formation of molecular complexes (dimers)
- Dissociation of complexes
- Combined translational and rotational diffusion

The simulation operates within a defined box with customizable physical parameters, making it suitable for studying a wide range of biological phenomena at the single-molecule level.

## Simulation Model

The diffusion simulation is based on the Smoluchowski dynamics model with the following components:

1. **Particles**: Represented as point particles (monomers) or rigid structures (dimers)
2. **Diffusion**: Isotropic Brownian motion with specified diffusion coefficients
3. **Reactions**: 
   - Association: Two monomers within reaction radius form a dimer
   - Dissociation: Dimers break with rate k_off
4. **Boundaries**: Periodic or reflecting boundary conditions

At each time step, the simulation:
- Updates molecular states (dimerization/dissociation)
- Updates positions with appropriate diffusion models
- Handles boundary conditions

### Physical Units

All simulation parameters use consistent physical units:
- Spatial dimensions: microns (μm)
- Time: seconds (s)
- Diffusion coefficients: μm²/s
- Rate constants: s⁻¹

## Getting Started

### Running a Basic Simulation

The main interface for running diffusion simulations is the `simulate` function with `DiffusionSMLMConfig`:

```julia
using SMLMSim

# Set simulation parameters
params = DiffusionSMLMConfig(
    density = 0.5,        # molecules per μm²
    box_size = 10.0,      # μm
    diff_monomer = 0.1,   # μm²/s
    diff_dimer = 0.05,    # μm²/s
    k_off = 0.2,          # s⁻¹
    r_react = 0.01,       # μm
    d_dimer = 0.05,       # μm
    dt = 0.01,            # s
    t_max = 10.0          # s
)

# Run the simulation
smld, info = simulate(params)
```

The `smld` output is a `BasicSMLD` structure containing all emitters across all time points, with each emitter having frame information corresponding to the camera settings. The `info` struct contains additional simulation metadata such as `info.elapsed_s`.

## Simulation Parameters

The `DiffusionSMLMConfig` structure allows you to customize various aspects of the simulation:

```julia
# More complex simulation
params = DiffusionSMLMConfig(
    density = 2.0,           # Higher density
    box_size = 20.0,         # Larger area
    diff_monomer = 0.2,      # Faster monomer diffusion
    diff_dimer = 0.08,       # Faster dimer diffusion
    diff_dimer_rot = 1.0,    # Faster rotational diffusion
    k_off = 0.05,            # Slower dissociation (more stable dimers)
    r_react = 0.015,         # Larger reaction radius
    d_dimer = 0.08,          # Larger dimer separation
    dt = 0.005,              # Smaller time step (higher precision)
    t_max = 30.0,            # Longer simulation
    ndims = 2,               # 2D simulation (default)
    boundary = "periodic",   # Periodic boundaries (default)
    camera_framerate = 20.0, # Camera frames per second
    camera_exposure = 0.04   # Camera exposure time
)
```

`dt` is the physics step and also sets the sub-steps per frame (motion blur):
`camera_exposure` and `1/camera_framerate` should be integer multiples of `dt`
(otherwise `simulate` rounds to the nearest step count and warns).
The `γ` argument of `simulate` is the emission rate in photons/s; each of the
`n_sub` records in a frame carries `γ·dt`, so a frame holds `γ·n_sub·dt` photons. That
equals `γ·camera_exposure` when `camera_exposure` is an integer multiple of `dt` and not longer
than the frame period `1/camera_framerate`; otherwise `n_sub = round(camera_exposure/dt)`
(at least 1), capped at the number of steps in the frame period, with a warning. Without `γ`, each record carries 1000 photons (γ = 1000/dt,
0.7's default; 0.8.0 will change the default to a fixed rate). The `photons` keyword is
deprecated (γ = photons/dt) and is removed in 0.8.0.

To continue a run, pass `starting_conditions=smld` or `extract_end_state(smld)`; both resume
at the exact end state of an unchanged run. Brightness: a run whose rate was set with `γ` continues at
that rate (each record carries γ·dt at the new `dt`) when every resumed molecule still carries γ·dt at
the saved `dt`, to a relative 1e-6, so a frame's brightness is unchanged as long as the effective
exposure `n_sub·dt` is unchanged; a new `dt` that changes `n_sub·dt` (an exposure that is not a whole
number of steps, or one capped at the frame period) changes the brightness with it. Otherwise (a
default or `photons` source, an SMLD with edited photons, or one of unknown provenance) each molecule
keeps its photons per record, as in 0.7.1, with a warning when a γ rate could not be kept. D: a track
keeps its saved D and mobility class (drawn from a `monomer_mobility` mixture) when the run saved one
for it, also in a filtered, time-cut or edited subset, with a warning when the new mixture differs from
the run's; every other track draws from the new mixture or uses `diff_monomer` at run time, which is
never saved, so a later change applies (such a track has no entry in `metadata["monomer_D"]`).
Continuation assumes the SMLD comes from one simulation run, or a filtered subset of one; continuing a
concatenation of different runs is unsupported, and the γ check cannot detect a molecule from another
run that happens to carry γ·dt. A concatenation or merge (SMLMData's `cat_smld` and `merge_smld` mark
it in the metadata) or, as a backstop, a last frame holding two records of one track at the same
timestamp, which no single run produces, is taken as unknown provenance: no γ, rate source or saved D
is carried, with one warning. One limitation: an SMLD that was filtered or edited resumes from each track's latest record in
its last frame, which is not the exact end state and carries no per-molecule history beyond that record
(blinking, bleaching or brightness-jitter state is not rebuilt). `extract_final_state` is deprecated.
The full continuation rule is in [Placement and Continuation Rules](rules.md).

Monomers can be given a mixture of mobility populations with
`monomer_mobility = [(0.85, 0.0), (0.05, 0.08), (0.10, 0.38)]` (entries are
`(fraction, D)`); the drawn coefficient per molecule is stored in
`smld.metadata["monomer_D"]`. `frame_dimer_truth(smld)` returns per-frame, per-molecule
dimer ground truth; its `mixed` field is `true` when exactly one of a molecule and its partner is
immobile (monomer D = 0).

A bound pair diffuses with `diff_dimer` by default (`pair_mobility = :fixed`, the 0.7 behaviour).
With `pair_mobility = :min` a pair moves at `min(D1, D2) × diff_dimer/diff_monomer`, where `D1`, `D2`
are its partners' monomer D, so two partners at `diff_monomer` move at `diff_dimer` (as under `:fixed`)
and a pair with an immobile partner (D = 0) does not move or rotate while bound: the immobile partner
keeps its position when the pair forms, and the mobile partner is placed `d_dimer` from it (in a
reflecting box, mirrored across the immobile partner on each axis it would leave, so the bond keeps its
length; if some axis fits neither way, possible only when `box_size < 2·d_dimer`, it keeps its
position). With `diff_dimer = 0` a `:min` pair does not rotate, while a `:fixed` pair still rotates
with `diff_dimer_rot`.

In a reflecting box a mobile pair (either setting) is a rigid body: when an end would leave the box at
formation, both partners move to `d_dimer/2` either side of their midpoint, shifted inward just enough
to fit, and while bound the pair's center folds off the walls moved in by each end's half-extent as
often as it crosses them, so both partners stay inside the box at `d_dimer` apart. When
`box_size < d_dimer` each partner is reflected on its own, as in 0.7.1. Under periodic boundaries (the
default) a bound pair moves from its partner's minimum image, so a pair straddling the boundary moves by
one step, and a pair forms from its partner's minimum image too (two monomers within `r_react` across the edge form a pair, which 0.7.1 never did), with each partner wrapped into the box on the formation step. Under the default `:fixed`, forming a pair places both partners `d_dimer` apart about their
midpoint and a bound pair moves with `diff_dimer` and rotates with `diff_dimer_rot`, even when a member
is immobile (a D = 0 population); `simulate` warns once per run when a step moves an immobile member,
recorded by the camera or not, and `pair_mobility = :min` keeps such a pair in place.

On dissociation the partners are placed at least `r_react` apart along the pair axis (the minimum-image distance under periodic boundaries), so a pair does not re-form at the next step only because `d_dimer < r_react`; a mobile partner can still diffuse back and re-form (geminate re-encounter). Under `pair_mobility = :min` an immobile partner stays put and the other is placed from it, mirrored across it on each axis it would leave a reflecting box. Otherwise, including every pair under `:fixed`, both partners move apart about their midpoint, shifted inward just enough to fit a reflecting box or wrapped under periodic boundaries; under `:fixed` this moves an immobile partner, and `simulate` counts it toward its one warning per run. Partners already `r_react` or more apart (by minimum image, computed in Float64 from their stored coordinates) keep their positions. That is every pair bound at `d_dimer` when `d_dimer` exceeds `r_react` by more than coordinate rounding and the box is at least `2·d_dimer`; in a smaller box a pair can be closer than `d_dimer` (a reflecting anchored formation that kept its separation, or a periodic bond whose minimum image is shorter), and one closer than `r_react` moves at the split.
The one exception to "an immobile partner stays put" under `:min`: when both partners are immobile, the one with the higher `track_id` is moved, because a pair closer than `r_react` must be separated to `r_react` and neither partner could otherwise ever move apart. When the placement cannot fit along the pair's axis (placed from an immobile partner, some axis fits neither way, possible only when `box_size` is under twice the separation; about their midpoint, some component of the separation along the axis exceeds `box_size`, or `box_size/2` under periodic boundaries), the partners stay where they are (the 0.7.1 behaviour), and `simulate` warns once per run. The separation is `r_react` plus a margin that survives rounding to the coordinate type. Each clause of pair placement, dissociation and continuation is stated in [Placement and Continuation Rules](rules.md).


## Microscope Image Generation

The diffusion simulation can be converted into realistic microscope images using point spread function models.

### Creating Images from Simulation Data

```julia
# Set up camera and PSF
pixelsize = 0.1  # 100nm pixels
pixels = Int64(round(params.box_size/pixelsize))
camera = IdealCamera(1:pixels, 1:pixels, pixelsize)

# Set up PSF (Gaussian with 150nm width)
using MicroscopePSFs
psf = MicroscopePSFs.GaussianPSF(0.15)  # 150nm PSF width

# Generate images
image_stack, img_info = gen_images(smld, psf;
    photons=1000.0,
    bg=5.0,
    poisson_noise=true
)
```

## Analyzing Results

### Extracting Dimers

To focus on the behavior of dimers, you can extract only the molecules in dimer state:

```julia
# Extract dimers from the full simulation
dimer_smld = get_dimers(smld)
```

### Analyzing Dimer Formation

To quantify the formation of dimers over time:

```julia
# Calculate fraction of molecules in dimer state
frames, dimer_fractions = analyze_dimer_fraction(smld)

# Plot dimer formation over time
using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1],
    xlabel="Frame",
    ylabel="Fraction of molecules in dimers",
    title="Dimer formation dynamics"
)
lines!(ax, frames, dimer_fractions)
fig
```

## Emitter Types

The diffusion module introduces specialized emitter types:

- `DiffusingEmitter2D`: 2D emitter with state information (monomer/dimer)
- `DiffusingEmitter3D`: 3D emitter with state information (monomer/dimer)

These types include additional properties:
- `timestamp`: Actual simulation time
- `state`: Molecular state (`:monomer` or `:dimer`)
- `partner_id`: ID of linked molecule (for dimers)

