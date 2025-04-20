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

## Kinetic Models

### Continuous Time Markov Chain

Fluorophore state transitions are modeled using a Continuous Time Markov Chain (CTMC):

```julia
# Create a CTMC for a two-state system
# State 1: ON (fluorescent), State 2: OFF (dark)
q = [-5 5; 10 -10]  # Units: s⁻¹
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
fluor = GenericFluor(γ=10000.0, q=[-5 5; 10 -10])
photons = intensity_trace(fluor, 1000, 50.0)
```

The `intensity_trace` function:
1. Simulates the CTMC state trajectory
2. Integrates photon emission (rate γ) during ON states
3. Accumulates photons within each frame's exposure time

## Common Photophysical Models

### Two-State Model

The simplest model contains just ON and OFF states:

```julia
# Two-state model (ON ⟷ OFF)
# kon = 5 s⁻¹, koff = 10 s⁻¹
fluor = GenericFluor(
    γ=1e4,                # 10,000 photons/s
    q=[-10 10; 5 -5]      # [ON→OFF; OFF→ON] rates in s⁻¹
)
```

This produces exponentially distributed ON and OFF times.

### Three-State Model with Bleaching

For more realistic behavior including photobleaching:

```julia
# Three-state model (ON ⟷ OFF → BLEACHED)
# State 1: ON, State 2: OFF, State 3: BLEACHED
fluor = GenericFluor(
    γ=1e4,
    q=[-10.1 10 0.1; 5 -5 0; 0 0 0]  # Note: state 3 is absorbing (no outgoing transitions)
)
```

The third state is irreversible (absorbing state), representing photobleaching.

## Using Photophysics in Simulations

The `simulate()` function integrates these photophysical models automatically:

```julia
# Simulation with custom fluorophore
camera = IdealCamera(128, 128, 0.1)
fluor = GenericFluor(γ=2e4, q=[-20 20; 5 -5])

smld_true, smld_model, smld_noisy = simulate(
    molecule=fluor,
    framerate=50.0,     # frames per second
    nframes=2000,       # total frames
    minphotons=100,     # detection threshold
    camera=camera
)
```

Behind the scenes, this uses the `kinetic_model()` function to apply the photophysical model to each emitter position.

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

