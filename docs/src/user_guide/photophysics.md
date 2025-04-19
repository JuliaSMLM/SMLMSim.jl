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
# Create a fluorophore with default parameters
fluor = GenericFluor()

# Create a fluorophore with custom parameters
fluor = GenericFluor(
    γ=1e5,                # photon emission rate in Hz
    q=[0 10; 1e-2 0]      # rate matrix in s⁻¹
)
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
q = [0 5; 10 0]  # Units: s⁻¹
simulation_time = 10.0  # seconds
initial_state = 2  # Start in dark state

ctmc = CTMC(q, simulation_time, initial_state)
```

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
fluor = GenericFluor(γ=10000.0, q=[0 5; 10 0])
photons = intensity_trace(fluor, 1000, 50.0)

# Plot the intensity trace
using CairoMakie
fig = Figure(size=(700, 300))
ax = Axis(fig[1, 1], xlabel="Frame", ylabel="Photons", title="Intensity Trace")
lines!(ax, 1:length(photons), photons)
fig
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
    q=[0 10; 5 0]         # [ON→OFF; OFF→ON] rates in s⁻¹
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
    q=[0 10 0.1; 5 0 0; 0 0 0]  # Note: state 3 is absorbing (no outgoing transitions)
)
```

The third state is irreversible (absorbing state), representing photobleaching.

## Using Photophysics in Simulations

The `simulate()` function integrates these photophysical models automatically:

```julia
# Simulation with custom fluorophore
camera = IdealCamera(1:128, 1:128, 0.1)
fluor = GenericFluor(γ=2e4, q=[0 20; 5 0])

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

## Customizing Photon Counts

The number of detected photons depends on several factors:

1. **Emission rate (γ)**: Base photon emission rate of the fluorophore
2. **Exposure time**: Determined by the framerate (1/framerate in seconds)
3. **State occupancy**: Time spent in the ON state during exposure
4. **Detection threshold**: Minimum photons for detection (minphotons)

You can adjust these parameters to simulate different experimental conditions:

```julia
# Bright fluorophore with high photon counts
bright_fluor = GenericFluor(γ=5e4, q=[0 5; 1 0])

# Dim fluorophore with low photon counts
dim_fluor = GenericFluor(γ=5e3, q=[0 10; 2 0])

# Simulate with different detection thresholds
smld_true, smld_model, smld_noisy = simulate(
    molecule=bright_fluor,
    minphotons=200,  # Higher threshold for good signal-to-noise
    framerate=100.0  # Faster frame rate = less time to collect photons
)
```

## Advanced: Multi-State Models

For more complex photophysical behavior, you can create larger rate matrices:

```julia
# Four-state model with multiple dark states
# State 1: ON, States 2-3: OFF (different lifetimes), State 4: BLEACHED
fluor = GenericFluor(
    γ=1e4,
    q=[0 5 1 0.01;   # ON -> OFF1, OFF2, BLEACHED
       2 0 0.5 0;    # OFF1 -> ON, OFF2
       0.2 0.1 0 0;  # OFF2 -> ON, OFF1
       0 0 0 0]      # BLEACHED (absorbing state)
)
```

This approach can model complex photophysics like triplet states, dark states with different lifetimes, and other fluorophore-specific behaviors.