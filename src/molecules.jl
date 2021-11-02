```
    Molecule

Photophysical properties of a molecule. 


This is the most general type of luminecent or scattering single molecule.  
Inherited types will defines the properties of a class of molecules. 

```
abstract type Molecule end 



```
    GenericFluor

Defines a fluorophore

# Fields
- `γ`: photon emission rate in Hz
- `q': state transision matrix. Default: q=[1.0]
```
mutable struct GenericFluor <: Molecule 
    γ::AbstractFloat         
    q::Array{AbstractFloat}=[1.0]
end

