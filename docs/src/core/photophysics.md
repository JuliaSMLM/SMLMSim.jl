```@meta
CurrentModule = SMLMSim
```

# Fluorophore Photophysics

This page covers the photophysical models in SMLMSim and how they're used to generate realistic blinking behavior for SMLM simulations.

## Overview

Fluorophores in single-molecule microscopy exhibit complex photophysical behavior, including:

- Stochastic switching between fluorescent (ON) and non-fluorescent (OFF) states
- Varying photon emission rates
- Photobleaching (permanent transition to a dark state)

SMLMSim models these behaviors using kinetic state models and stochastic simulation of state transitions.

## Fluorophore Models

### GenericFluor

The main fluorophore model in SMLMSim is the `GenericFluor` type:

```julia
using SMLMSim

# Define a two-state fluorophore using the positional constructor
# Note: q matrix already had correct diagonal elements, just removed γ=
fluor = GenericFluor(10000.0, [-5.0 5.0; 10.0 -10.0])

# ... rest of the example ...
```

Parameters:
- `γ`: Photon emission rate in Hz (photons per second)
- `q`: State transition rate matrix where `q[i,j]` is the transition rate from state i to j (in s⁻¹)

By convention, state 1 is the fluorescent (ON) state, and other states are non-fluorescent (OFF).

The keyword constructor accepts real numbers, including integers, for `photons`, `k_off`,
and `k_on`. Floating point values retain their types and precision;
only non-float reals such as integers and rationals are converted to Float64:

```julia
fluor = GenericFluor(; photons=10000, k_off=50, k_on=1)
```

Through 0.7.3 the keyword values had to be `AbstractFloat`, so these integer values threw a
type error. The positional constructor preserves the types of the supplied rate matrix
and emission rate.

## Kinetic Models

### Continuous Time Markov Chain

Fluorophore state transitions are modeled using a Continuous Time Markov Chain (CTMC):

```julia
# Create a CTMC for a two-state system
# State 1: ON (fluorescent), State 2: OFF (dark)
q = [-5.0 5.0; 10.0 -10.0]  # Units: s⁻¹
simulation_time = 10.0  # seconds
initial_state = 2  # Start in dark state

ctmc = CTMC(q, simulation_time, initial_state)
```

#### Understanding the Rate Matrix

The transition rate matrix `q` represents the rates at which the system transitions between states:

- `q[i,j]` (i≠j): Rate of transition from state i to state j
- `q[i,i]`: Negative sum of all outgoing rates from state i

For example, in a two-state system:

```
q = [-k_off  k_off;
      k_on   -k_on]
```

Where:
- `k_off`: Rate of transitioning from ON (state 1) to OFF (state 2)
- `k_on`: Rate of transitioning from OFF (state 2) to ON (state 1)
- Diagonal elements are negative sums of their respective rows

The CTMC simulates transitions between states by:
1. Sampling the time until the next transition (exponentially distributed with rate -q[current,current])
2. Sampling the next state with probabilities proportional to the transition rates

The CTMC provides a complete trajectory of state transitions:

```julia
# Get state at a specific time
state_at_1s = get_state(ctmc, 1.0)

# Get next state transition after a specific time
next_state, transition_time = get_next(ctmc, 0.5)
```

### Intensity Traces

To simulate the fluorescence signal over time, SMLMSim integrates photon emission during ON states:

```julia
# Generate intensity trace for 1000 frames at 50 fps
fluor = GenericFluor(10000.0, [-5.0 5.0; 10.0 -10.0])
photons = intensity_trace(fluor, 1000, 50.0)
```

The `intensity_trace` function:
1. Simulates the CTMC state trajectory
2. Integrates photon emission (rate γ) during ON states
3. Accumulates photons within each frame's exposure time

When building the chain, the rate matrix and simulation time are promoted to a common
floating point type. For example, `GenericFluor(1f4, Float32[-50 50; 1f-2 -1f-2])` works
with the usual Float64 time arguments. Through 0.7.3 this combination failed with a method
error. Float64 simulations keep their existing random draws and photon counts.

An absorbing state holds through the frame end. Through 0.7.3 a bleached molecule could hang
when rounding placed the final frame end just past the CTMC's end time.

## Common Photophysical Models

### Two-State Model

The simplest model contains just ON and OFF states:

```julia
# Two-state model (ON ⟷ OFF)
# kon = 5 s⁻¹, koff = 10 s⁻¹
fluor = GenericFluor(
    1e4,                  # 10,000 photons/s
    [-10.0 10.0; 5.0 -5.0]  # [ON→OFF; OFF→ON] rates in s⁻¹
)
```

This produces exponentially distributed ON and OFF times.

### Three-State Model with Bleaching

For more realistic behavior including photobleaching:

```julia
# Three-state model (ON ⟷ OFF, ON → BLEACHED)
# State 1: ON, State 2: OFF, State 3: BLEACHED
fluor = GenericFluor(
    1e4,
    [-10.1 10.0 0.1; 5.0 -5.0 0.0; 0.0 0.0 0.0]    # Note: state 3 is absorbing (no outgoing transitions)
)
```

The third state is irreversible (absorbing state), representing photobleaching.

## Using Photophysics in Simulations

The `simulate()` function integrates these photophysical models automatically:

```julia
# Simulation with custom fluorophore
camera = IdealCamera(128, 128, 0.1)
fluor = GenericFluor(2e4, [-20.0 20.0; 5.0 -5.0])

params = StaticSMLMConfig(
    framerate=50.0,     # frames per second
    nframes=2000,       # total frames
    minphotons=100      # detection threshold
)
smld_noisy, info = simulate(params; molecule=fluor, camera=camera)
```

Behind the scenes, this uses the `kinetic_model()` function to apply the photophysical model to each emitter position.

`kinetic_model` uses the first record for each non-zero `track_id`, in sorted id order, and
warns once if later records for that id have different positions. Each zero id is a separate
unlinked molecule, processed afterward in input order and assigned an output id above the
largest non-zero input id. Through 0.7.3 all zero ids were grouped as one molecule, and
conflicting linked positions were dropped silently.

Input containers may have an abstract element type or mix numeric precisions within one
dimension. The 2D or 3D output type is chosen from the individual fits and its coordinates
are Float64. Mixing 2D and 3D inputs throws `ArgumentError`.

## Molecule-to-Molecule Brightness

By default every molecule emits at the model's rate `γ`, so every full-frame molecule carries
exactly `γ / framerate` photons. The keyword `brightness` of `simulate` (static) and
`kinetic_model` takes any univariate distribution of a dimensionless factor `b`. A molecule
with factor `b` emits at `b * γ` for all of its frames; the rate matrix is unchanged. One factor
is drawn per molecule and dataset, and the factor must not be negative. Redrawing per dataset is
intended: each dataset is an independent sample, as its initial states and blinking already are,
so the same position can differ in brightness between datasets.

For a mean and standard deviation of photons per full frame with a floor, set `γ` from the mean
and give the factor as a normal of mean 1, kept at or above the floor:

```julia
using SMLMSim, Distributions   # Normal, truncated and censored come from Distributions.jl

mean_ph, sd_ph, floor_ph = 3000.0, 600.0, 1000.0   # photons per full frame
fluor = GenericFluor(mean_ph * 50.0, [-10.0 10.0; 0.5 -0.5])   # γ for a 50 Hz frame rate
params = StaticSMLMConfig(framerate=50.0)

# cut at the floor and renormalized: no draw lies at the floor
b_trunc = truncated(Normal(1, sd_ph / mean_ph); lower=floor_ph / mean_ph)
# draws below the floor are set to it, as DECODE does
b_cens = censored(Normal(1, sd_ph / mean_ph); lower=floor_ph / mean_ph)

smld, info = simulate(params; molecule=fluor, brightness=b_trunc)
info.smld_model.metadata["brightness_factors"]   # (dataset, track_id) => factor
```

The two differ in what they do with the probability below the floor. `truncated` removes it
and renormalizes, so the distribution is smooth and no molecule sits exactly at the floor.
`censored` moves it onto the floor, so a fraction `cdf(Normal(1, sd/mean), floor/mean)` of the
molecules has the factor exactly equal to the floor. Both raise the mean above that of the
underlying normal and narrow its spread; `mean` and `var` of either distribution give the exact
moments. With 32-bit parameters, `truncated(Normal(1f0, 0.3f0); lower=0.7f0)` fails when drawn
from in Distributions 0.25; build the distribution from 64-bit values.

A distribution with negative values (such as a plain `Normal`) throws `ArgumentError`. With
`brightness` given, `metadata["brightness"]` holds the distribution and
`metadata["brightness_factors"]` the factor drawn for each `(dataset, track_id)`; without it no
random number is drawn for brightness and neither key is present.

## Photophysical Parameters and Duty Cycle

The photophysical properties determine the "duty cycle" - the fraction of time a fluorophore is in the ON state:

For a two-state model, the duty cycle is:
```
duty_cycle = kon / (kon + koff)
```

Where:
- `kon` is the OFF→ON rate (q[2,1] in the rate matrix)
- `koff` is the ON→OFF rate (q[1,2] in the rate matrix)

Typical duty cycles for SMLM fluorophores range from 0.0001 to 0.01.

