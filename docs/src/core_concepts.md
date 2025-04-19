```@meta
CurrentModule = SMLMSim
```

# Core Concepts

This page introduces the fundamental concepts behind Single Molecule Localization Microscopy (SMLM) simulation and how SMLMSim implements them.

## Physical Units

All simulations in SMLMSim use physical units:
- Spatial coordinates: microns (μm)
- Time: seconds (s)
- Rates: per second (s⁻¹)

This consistent use of physical units allows direct comparison with experimental data and facilitates integration with other physical models.

## Coordinate Systems

- **2D simulations**: Planar (x,y) coordinates, typically representing the focal plane
- **3D simulations**: Volume (x,y,z) coordinates, with z representing axial position
- **Origin**: (0,0) or (0,0,0) typically represents the center of the field of view
- **Camera pixels**: Defined by physical size (e.g., 100nm = 0.1μm)

## Molecular Patterns

Patterns define the spatial arrangement of fluorophores:

- **Basic patterns**: 
  - `Nmer2D`/`Nmer3D`: N molecules arranged in a circle (common for protein complexes)
  - `Line2D`/`Line3D`: Linear arrangements with specified density

- **Pattern properties**:
  - Diameter/length: Physical size in microns
  - Position: Coordinates within the field of view
  - Orientation: Rotational state (can be randomized)

- **Pattern distribution**:
  - Density (ρ): Number of patterns per square micron
  - Random positioning within the field of view
  - Random rotational orientation

## Fluorophore Photophysics

Fluorophores in SMLMSim are modeled using stochastic kinetic models:

- **States**: Typically ON (fluorescent) and OFF (dark) states
- **Transitions**: Characterized by rate constants (in s⁻¹)
- **Photon emission**: Occurs during ON states with specified rate γ (in Hz)
- **Continuous Time Markov Chain (CTMC)**: Used to simulate state transitions

The photophysical behavior is used for generating realistic blinking patterns observed in SMLM experiments.

## Simulation Pipeline

A typical simulation follows these steps:

1. **Pattern generation**: Create spatial distribution of fluorophores
2. **Kinetic modeling**: Apply stochastic blinking behavior
3. **Localization uncertainty**: Add noise based on photon statistics
4. **Image formation** (optional): Generate microscope images using PSF models

The high-level `simulate()` function integrates these steps into a convenient interface.

## Output Data Structure

Simulation results are organized in `SMLD` (Single Molecule Localization Data) structures:

- **Ground truth**: True molecular positions without blinking or uncertainty
- **Kinetic model**: Subset of true positions that appear in each frame based on blinking
- **Noisy data**: Realistic SMLM data with both blinking and position errors

These structures contain individual emitters with properties such as position, frame number, photon count, and uncertainty estimates.

## Interaction-Diffusion Models

Beyond static patterns, SMLMSim can model dynamic processes:

- **Diffusion**: Brownian motion of molecules with specified diffusion coefficients
- **Interaction**: Formation and dissociation of molecular complexes
- **Smoluchowski dynamics**: Physically accurate simulation of reaction-diffusion processes

These advanced capabilities allow simulation of dynamic biological processes and single-particle tracking experiments.