"""
    Molecule

Photophysical properties of a molecule. 


This is the most general type of luminecent or scattering single molecule.  
Inherited types will defines the properties of a class of molecules. 

"""
abstract type Molecule end



"""
    GenericFluor

Defines a fluorophore

# Fields
- γ: photon emission rate in Hz, Default: 1e3
- q: state transision matrix. Default: q=[1.0]
"""
mutable struct GenericFluor <: Molecule
    γ::AbstractFloat
    q::Array{<:AbstractFloat}
end
function GenericFluor(;
    γ::AbstractFloat=1e5,
    q::Array{<:AbstractFloat}=[1]
)
    return GenericFluor(γ, q)
end








