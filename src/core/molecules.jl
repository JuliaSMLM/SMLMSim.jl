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

Positional construction preserves the floating point types. `intensity_trace` promotes the
rate matrix and simulation time to a common floating point type when building the chain.
Through 0.7.3 a Float32 rate matrix with a Float64 time failed with a method error.

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
```
"""
struct GenericFluor <: Molecule
    γ::AbstractFloat
    q::Array{<:AbstractFloat}
end


"""
    GenericFluor(; photons::Real=1e5, k_off::Real=50.0, k_on::Real=1e-2)

Create a simple two-state (on/off) fluorophore with specified parameters.

# Arguments
- `photons::Real`: Photon emission rate in photons/s
- `k_off::Real`: Off-switching rate (on→off) in s⁻¹
- `k_on::Real`: On-switching rate (off→on) in s⁻¹

# Details
Creates a fluorophore with a 2-state model and the specified rates.
State 1 is the on (bright) state, and state 2 is the off (dark) state.
The rate matrix is constructed as: q = [-k_off k_off; k_on -k_on]
`AbstractFloat` keyword values retain their types and precision. Only
non-float reals (such as integers and rationals) are converted to Float64; the rate matrix
uses normal numeric promotion. Through 0.7.3 these keywords required `AbstractFloat`
values, so valid integer values such as `photons=10000` threw a type error.

Note: k_on and k_off are transition rates (1/s), not duty cycle fractions.
The duty cycle (fraction of time in ON state) is k_on/(k_on + k_off).
For typical dSTORM, k_on << k_off gives low duty cycle (mostly dark, brief blinks).
"""
function GenericFluor(;
    photons::Real=1e5,
    k_off::Real=50.0,
    k_on::Real=1e-2
)
    photons = photons isa AbstractFloat ? photons : Float64(photons)
    k_off = k_off isa AbstractFloat ? k_off : Float64(k_off)
    k_on = k_on isa AbstractFloat ? k_on : Float64(k_on)
    # Create rate matrix for 2-state system
    q = [-k_off k_off; k_on -k_on]
    
    return GenericFluor(photons, q)
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