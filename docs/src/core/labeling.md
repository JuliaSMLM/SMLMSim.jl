```@meta
CurrentModule = SMLMSim
```

# Labeling Strategies

Labeling controls how many fluorophores attach to each binding site in a pattern. This is separate from photophysics — labeling determines *how many* fluorophores are placed, while [`Molecule`](@ref) types handle the *blinking kinetics* of each one.

## Pipeline Position

Labeling sits between pattern generation and photophysics in the simulation pipeline:

```
Pattern (geometry) → uniform2D/3D() → binding site coordinates
  → apply_labeling() → fluorophore positions (one or more per site)
  → kinetic_model() → blinking behavior per fluorophore
  → apply_noise() → localization uncertainty
```

## Labeling Types

All labeling types support an `efficiency` parameter (0 to 1) controlling the probability that a binding site gets labeled at all. This is applied *before* sampling the number of fluorophores.

### FixedLabeling

Deterministic: exactly `n` fluorophores per labeled site.

```julia
# Default: 1 fluorophore per site, 100% efficiency
labeling = FixedLabeling()

# 2 fluorophores per site, 90% of sites labeled
labeling = FixedLabeling(2; efficiency=0.9)
```

### PoissonLabeling

Poisson-distributed number of fluorophores per site. Useful when the labeling count varies stochastically (e.g., random antibody binding).

```julia
# Average 1.5 fluorophores per site
labeling = PoissonLabeling(1.5)

# Average 2.0 per site, 80% labeling efficiency
labeling = PoissonLabeling(2.0; efficiency=0.8)
```

Note: with Poisson statistics, some sites may receive 0 fluorophores even when `efficiency=1.0`, especially for small mean values.

### BinomialLabeling

Binomial-distributed: each site has `n` potential attachment points, each occupied with probability `p`. Models scenarios like an antibody with multiple dye-conjugation sites.

```julia
# 4 attachment points, 80% occupancy each
labeling = BinomialLabeling(4, 0.8)

# Same, but only 90% of sites receive an antibody
labeling = BinomialLabeling(4, 0.8; efficiency=0.9)
```

## Using Labeling in Simulations

Pass a labeling strategy to `simulate()` via the `labeling` keyword:

```julia
params = StaticSMLMConfig(density=1.0, σ_psf=0.13)

# Default (FixedLabeling with 1 fluorophore per site)
smld_noisy, info = simulate(params; pattern=Nmer2D())

# Poisson labeling
smld_noisy, info = simulate(params; pattern=Nmer2D(), labeling=PoissonLabeling(1.5))

# Binomial labeling with reduced efficiency
smld_noisy, info = simulate(params;
    pattern=Nmer2D(n=6, d=0.2),
    labeling=BinomialLabeling(4, 0.8; efficiency=0.9)
)
```

## Direct Usage with apply_labeling

For custom workflows, use `apply_labeling` directly on binding site coordinates:

```julia
# Generate binding sites from a pattern distribution
x, y, pattern_ids = uniform2D(1.0, Nmer2D(), 10.0, 10.0)

# Apply labeling with pattern ID tracking
(x_labeled, y_labeled), new_ids = apply_labeling((x, y), pattern_ids, PoissonLabeling(1.5))

# Without pattern tracking (backward-compatible)
x_labeled, y_labeled = apply_labeling((x, y), FixedLabeling(2))

# Works with 3D coordinates
x, y, z, pattern_ids = uniform3D(1.0, Nmer3D(), 10.0, 10.0)
(x_labeled, y_labeled, z_labeled), new_ids = apply_labeling((x, y, z), pattern_ids, BinomialLabeling(4, 0.8))
```

## Sampling Fluorophore Counts Directly

Use `n_fluorophores` to sample the count for a single site:

```julia
labeling = PoissonLabeling(1.5)
n = n_fluorophores(labeling)  # Returns Int (can be 0)
```

See the [API Reference](@ref) page for complete docstrings.
