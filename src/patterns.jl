#define patterned stuctures

```
    Pattern

Abstract type for structured patterns of molecules    
```
abstract type Pattern end


```
    Nmer <: Pattern

```
mutable struct Nmer <: Pattern
    d::AbstractFloat
    
end


```
    Uniform <: Pattern

Generate data with uniform random distribution.    

# Fields
- `ρ`: molecule density

```
mutable struct Uniform <: Pattern
    ρ::AbstractFloat   
end




