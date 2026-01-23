# Labeling type definitions
"""
    AbstractLabeling

Abstract type for labeling strategies that determine how many fluorophores
attach to each binding site.

Labeling is separate from photophysics (Molecule) - this type hierarchy handles
the statistics of fluorophore attachment, while Molecule handles blinking kinetics.

# Implementations
- `FixedLabeling`: Deterministic number of fluorophores per site
- `PoissonLabeling`: Poisson-distributed number of fluorophores per site
- `BinomialLabeling`: Binomial-distributed number of fluorophores per site
"""
abstract type AbstractLabeling end

#==========================================================================
Concrete Labeling Types
==========================================================================#

"""
    FixedLabeling <: AbstractLabeling

Deterministic labeling with exactly `n` fluorophores per site.

# Fields
- `n::Int`: Number of fluorophores per site
- `efficiency::Float64`: Probability that a site gets labeled at all (0 to 1)

# Examples
```julia
# Perfect labeling - exactly 1 fluorophore per site (default/current behavior)
labeling = FixedLabeling()

# 2 fluorophores per site, 90% of sites labeled
labeling = FixedLabeling(2; efficiency=0.9)
```
"""
struct FixedLabeling <: AbstractLabeling
    n::Int
    efficiency::Float64

    function FixedLabeling(n::Int, efficiency::Float64)
        n < 0 && throw(ArgumentError("n must be non-negative"))
        (efficiency < 0 || efficiency > 1) && throw(ArgumentError("efficiency must be between 0 and 1"))
        new(n, efficiency)
    end
end

FixedLabeling(; n::Int=1, efficiency::Float64=1.0) = FixedLabeling(n, efficiency)
FixedLabeling(n::Int; efficiency::Float64=1.0) = FixedLabeling(n, efficiency)

function Base.show(io::IO, l::FixedLabeling)
    print(io, "FixedLabeling(n=$(l.n), efficiency=$(l.efficiency))")
end

"""
    PoissonLabeling <: AbstractLabeling

Poisson-distributed number of fluorophores per site.

# Fields
- `mean::Float64`: Mean number of fluorophores per site (λ for Poisson)
- `efficiency::Float64`: Probability that a site gets labeled at all (0 to 1)

# Examples
```julia
# Average 1.5 fluorophores per site
labeling = PoissonLabeling(1.5)

# Average 2 fluorophores per site, but only 80% of sites get labeled
labeling = PoissonLabeling(2.0; efficiency=0.8)
```

# Note
With Poisson statistics, some sites may receive 0 fluorophores even when
efficiency=1.0 (especially for small mean values). The efficiency parameter
controls a separate "does this site get labeled at all" step.
"""
struct PoissonLabeling <: AbstractLabeling
    mean::Float64
    efficiency::Float64

    function PoissonLabeling(mean::Float64, efficiency::Float64)
        mean < 0 && throw(ArgumentError("mean must be non-negative"))
        (efficiency < 0 || efficiency > 1) && throw(ArgumentError("efficiency must be between 0 and 1"))
        new(mean, efficiency)
    end
end

PoissonLabeling(; mean::Float64=1.0, efficiency::Float64=1.0) = PoissonLabeling(mean, efficiency)
PoissonLabeling(mean::Float64; efficiency::Float64=1.0) = PoissonLabeling(mean, efficiency)

function Base.show(io::IO, l::PoissonLabeling)
    print(io, "PoissonLabeling(mean=$(l.mean), efficiency=$(l.efficiency))")
end

"""
    BinomialLabeling <: AbstractLabeling

Binomial-distributed number of fluorophores per site. Models scenarios where
each binding site has `n` potential attachment points, each occupied with
probability `p`.

# Fields
- `n::Int`: Number of potential attachment points per site
- `p::Float64`: Probability each attachment point is occupied (0 to 1)
- `efficiency::Float64`: Probability that a site gets labeled at all (0 to 1)

# Examples
```julia
# Antibody with 4 dye attachment points, each 80% likely to be occupied
labeling = BinomialLabeling(4, 0.8)

# Same, but only 90% of binding sites get an antibody
labeling = BinomialLabeling(4, 0.8; efficiency=0.9)
```
"""
struct BinomialLabeling <: AbstractLabeling
    n::Int
    p::Float64
    efficiency::Float64

    function BinomialLabeling(n::Int, p::Float64, efficiency::Float64)
        n < 0 && throw(ArgumentError("n must be non-negative"))
        (p < 0 || p > 1) && throw(ArgumentError("p must be between 0 and 1"))
        (efficiency < 0 || efficiency > 1) && throw(ArgumentError("efficiency must be between 0 and 1"))
        new(n, p, efficiency)
    end
end

BinomialLabeling(; n::Int=1, p::Float64=1.0, efficiency::Float64=1.0) = BinomialLabeling(n, p, efficiency)
BinomialLabeling(n::Int, p::Float64; efficiency::Float64=1.0) = BinomialLabeling(n, p, efficiency)

function Base.show(io::IO, l::BinomialLabeling)
    print(io, "BinomialLabeling(n=$(l.n), p=$(l.p), efficiency=$(l.efficiency))")
end

#==========================================================================
Interface Functions
==========================================================================#

"""
    n_fluorophores(labeling::AbstractLabeling) -> Int

Sample the number of fluorophores to place at a single binding site.

This function is called once per binding site during labeling. It first
checks labeling efficiency (probability site gets labeled at all), then
samples from the appropriate distribution.

# Returns
- `Int`: Number of fluorophores (can be 0)
"""
function n_fluorophores(l::FixedLabeling)
    rand() > l.efficiency && return 0
    return l.n
end

function n_fluorophores(l::PoissonLabeling)
    rand() > l.efficiency && return 0
    return rand(Poisson(l.mean))
end

function n_fluorophores(l::BinomialLabeling)
    rand() > l.efficiency && return 0
    return rand(Binomial(l.n, l.p))
end

#==========================================================================
Apply Labeling to Coordinates
==========================================================================#

"""
    apply_labeling(coords, labeling::AbstractLabeling)

Apply labeling strategy to binding site coordinates, expanding each site
to the appropriate number of fluorophore positions.

# Arguments
- `coords`: Tuple of coordinate vectors `(x, y)` for 2D or `(x, y, z)` for 3D
- `labeling::AbstractLabeling`: Labeling strategy to apply

# Returns
- Tuple of coordinate vectors with fluorophore positions

# Example
```julia
# Original: 100 binding sites
x, y = uniform2D(1.0, Nmer2D(), 10.0, 10.0)

# After Poisson labeling: variable number of fluorophores
x_labeled, y_labeled = apply_labeling((x, y), PoissonLabeling(1.5))
```
"""
function apply_labeling(coords::Tuple{Vector{T}, Vector{T}}, labeling::AbstractLabeling) where T
    x, y = coords
    new_x = T[]
    new_y = T[]

    for i in eachindex(x)
        n = n_fluorophores(labeling)
        for _ in 1:n
            push!(new_x, x[i])
            push!(new_y, y[i])
        end
    end

    return (new_x, new_y)
end

function apply_labeling(coords::Tuple{Vector{T}, Vector{T}, Vector{T}}, labeling::AbstractLabeling) where T
    x, y, z = coords
    new_x = T[]
    new_y = T[]
    new_z = T[]

    for i in eachindex(x)
        n = n_fluorophores(labeling)
        for _ in 1:n
            push!(new_x, x[i])
            push!(new_y, y[i])
            push!(new_z, z[i])
        end
    end

    return (new_x, new_y, new_z)
end
