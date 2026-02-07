# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build, Test & Run Commands

```bash
# Install dependencies
julia --project -e "using Pkg; Pkg.instantiate()"

# Run all tests
julia --project -e "using Pkg; Pkg.test()"

# Run specific test file (matches testset names like "Patterns", "Static", "Core")
julia --project -e "using Pkg; Pkg.test(\"SMLMSim\", test_args=[\"Patterns\"])"

# Build documentation
julia --project=docs/ docs/make.jl

# Run dev scripts
julia --project=dev/ dev/benchmark.jl
```

## Architecture Overview

SMLMSim simulates Single Molecule Localization Microscopy (SMLM) data. The package has four main modules that build on SMLMData.jl types:

```
SMLMSim (main module - re-exports from submodules)
├── Core           - Shared types and photophysics
├── StaticSMLM     - Fixed emitter simulation with blinking
├── InteractionDiffusion - Diffusing/interacting molecules (SPT)
└── CameraImages   - Generate camera images from emitters
```

### Simulation Pipeline (Static)

```
Pattern (geometry)
  → uniform2D/3D() → binding site coordinates
  → apply_labeling() → fluorophore positions (via AbstractLabeling)
  → create emitters → smld_true (BasicSMLD)
  → kinetic_model() → smld_model (with blinking, uses Molecule)
  → apply_noise() → smld_noisy (with localization uncertainty)
```

Key separation of concerns:
- **Pattern**: Spatial geometry of binding sites (Nmer2D, Line3D, etc.)
- **Labeling**: How many fluorophores per site (FixedLabeling, PoissonLabeling, BinomialLabeling)
- **Molecule**: Photophysics/blinking kinetics (GenericFluor with CTMC rate matrix)

### Type Hierarchy

```julia
# Patterns (geometry)
Pattern → Pattern2D → Nmer2D, Line2D
        → Pattern3D → Nmer3D, Line3D

# Labeling (fluorophore attachment statistics)
AbstractLabeling → FixedLabeling, PoissonLabeling, BinomialLabeling

# Molecules (photophysics)
Molecule → GenericFluor

# Simulation params
SMLMSimParams → StaticSMLMConfig, DiffusionSMLMConfig
```

### Key Files by Function

| Purpose | Files |
|---------|-------|
| Pattern geometry | `src/core/patterns.jl` |
| Labeling statistics | `src/core/labeling.jl` |
| Photophysics/blinking | `src/core/molecules.jl`, `src/core/photophysics.jl` |
| CTMC (state transitions) | `src/core/ctmc.jl` |
| Static simulation | `src/static/simulation.jl`, `src/static/parameters.jl` |
| Diffusion simulation | `src/diffusion/smoluchowski.jl`, `src/diffusion/types.jl` |
| Image generation | `src/camera_images/gen_images.jl`, `src/camera_images/noise.jl` |

### Dispatch Pattern

The `simulate()` function dispatches on parameter type:
- `simulate(params::StaticSMLMConfig; ...)` → static simulation (returns 3-tuple)
- `simulate(params::DiffusionSMLMConfig; ...)` → diffusion simulation (returns BasicSMLD)

## Code Style

- Physical units: μm for space, seconds for time
- Re-export SMLMData types rather than redefining
- Follow existing patterns when adding new Pattern, Labeling, or Molecule types
- Type parameters explicit for multiple dispatch functions
- Tests use `@test` with tolerances for floating-point comparisons
