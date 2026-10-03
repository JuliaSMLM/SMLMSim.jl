"""
    Molecule

Abstract type for representing photophysical properties of a molecule.

This is the most general type of luminescent or scattering single molecule.
Inherited types will define the properties of specific classes of molecules.
"""
abstract type Molecule end

"""
    GenericFluor <: Molecule

Defines a fluorophore with photophysical properties.

# Fields
- `γ::AbstractFloat`: Photon emission rate in Hz. Default: 1e5
- `q::Array{<:AbstractFloat}`: Rate matrix where q[i,j] for i≠j is the transition rate from state i to j,
  and q[i,i] is the negative exit rate from state i. Default: standard 2-state model with
  on->off rate of 50Hz and off->on rate of 1e-2Hz

# Examples
```julia
# Create a fluorophore with default parameters (using the 2-state keyword constructor)
fluor = GenericFluor()

# Create a fluorophore with custom parameters using the positional constructor
fluor = GenericFluor(1e5, [-50.0 50.0; 1e-2 -1e-2])

# Create a fluorophore using the 2-state keyword constructor
fluor = GenericFluor(; photons=1e5, k_off=10.0, k_on=1e-1)

# Create a fluorophore from a rate and a rate matrix by keyword
fluor = GenericFluor(; γ=1e4, q=[-10.0 10.0; 1e-1 -1e-1])
```
"""
struct GenericFluor <: Molecule
    γ::AbstractFloat
    q::Array{<:AbstractFloat}
end


"""
    GenericFluor(; photons, γ, k_off, k_on, q)

Create a fluorophore by keyword. Without `q` it is a simple two-state (on/off) fluorophore.

# Keywords
- `γ::Real` or `photons::Real`: Photon emission rate in Hz (default 1e5). Give at most one of the two.
- `q::AbstractMatrix`: Rate matrix, as in the positional constructor. Excludes `k_off` and `k_on`.
- `k_off::Real`: Off-switching rate (on→off) in Hz (default 50.0)
- `k_on::Real`: On-switching rate (off→on) in Hz (default 1e-2)

# Details
Without `q`, the 2-state rate matrix is `q = [-k_off k_off; k_on -k_on]`.
State 1 is the on (bright) state, and state 2 is the off (dark) state.

Note: k_on and k_off are transition rates (1/s), not duty cycle fractions.
The duty cycle (fraction of time in ON state) is k_on/(k_on + k_off).
For typical dSTORM, k_on << k_off gives low duty cycle (mostly dark, brief blinks).

Throws `ArgumentError` when both `γ` and `photons` are given, or when `q` is given with `k_off` or `k_on`.
"""
function GenericFluor(;
    photons::Union{Nothing,Real}=nothing,
    γ::Union{Nothing,Real}=nothing,
    k_off::Union{Nothing,Real}=nothing,
    k_on::Union{Nothing,Real}=nothing,
    q::Union{Nothing,AbstractMatrix{<:Real}}=nothing
)
    (γ === nothing || photons === nothing) ||
        throw(ArgumentError("give either γ or photons, not both"))
    (q === nothing || (k_off === nothing && k_on === nothing)) ||
        throw(ArgumentError("give either q or k_off/k_on, not both"))
    rate = float(γ !== nothing ? γ : photons !== nothing ? photons : 1e5)
    if q === nothing
        koff = float(k_off === nothing ? 50.0 : k_off)
        kon = float(k_on === nothing ? 1e-2 : k_on)
        q = [-koff koff; kon -kon]
    else
        q = Matrix(float.(q))
    end
    return GenericFluor(rate, q)
end

function Base.show(io::IO, fluor::GenericFluor)
    n_states = size(fluor.q, 1)
    print(io, "GenericFluor($(n_states) states, γ=$(fluor.γ) Hz)")
end

function Base.show(io::IO, ::MIME"text/plain", fluor::GenericFluor)
    n_states = size(fluor.q, 1)
    println(io, "GenericFluor with $n_states states:")
    println(io, "  Photon emission rate (γ) = $(fluor.γ) Hz")
    println(io, "  Rate matrix (q):")
    
    # Format the rate matrix with aligned columns
    for i in 1:n_states
        print(io, "    ")
        for j in 1:n_states
            val = fluor.q[i, j]
            # Format the rate value with appropriate precision
            val_str = abs(val) < 0.01 ? @sprintf("%.2e", val) : @sprintf("%.3f", val)
            print(io, lpad(val_str, 10))
        end
        if i < n_states
            println(io)
        end
    end
end