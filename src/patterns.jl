#define patterned stuctures

"""
    Pattern

Abstract type for structured patterns of molecules    
"""
abstract type Pattern end


"""
    Nmer2D <: Pattern

N molecules symmetricaly organized with diameter d    

Nmer2D(n::Int,d::AbstractFloat)

"""
mutable struct Nmer2D <: Pattern
    n::Int
    d::AbstractFloat
    x::Vector{AbstractFloat}
    y::Vector{AbstractFloat}
end
function Nmer2D(n::Int, d::AbstractFloat)

    nmer = Nmer2D(n, d, zeros(n), zeros(n))
    for nn = 1:n
        θ = 2 * pi / n * (nn - 1)
        nmer.x[nn] = d / 2 * cos(θ)
        nmer.y[nn] = d / 2 * sin(θ)
    end
    return nmer
end




"""
    Uniform <: Pattern

Generate data with uniform random distribution.    

# Fields
- `ρ`: molecule density

"""
mutable struct Point2D <: Pattern
    n::Int
    x::Vector{AbstractFloat}
    y::Vector{AbstractFloat}
    Point2D() = new(1, [0.0], [0.0])
end

"""
    function uniformPattern2D(ρ,xsize::AbstractFloat,ysize::AbstractFloat, p::Pattern)

create true positions of molecules from uniformly randomly placed patterns
"""
function uniform2D(ρ, p::Pattern, xsize::AbstractFloat, ysize::AbstractFloat)

    npatterns = rand(Poisson(xsize * ysize * ρ))
    ntotal = npatterns * p.n

    #make smd 
    smd = SMLMData.SMLD2D(ntotal)
    smd.datasize = Int.(ceil.([ysize; xsize]))
    for nn = 1:npatterns
        θ = 2 * pi * rand()
        x0=rand()*xsize
        y0=rand()*ysize
        
        for mm = 1:p.n
            idx = (p.n) * (nn - 1) + mm
            smd.x[idx] = p.x[mm] * cos(θ) - p.y[mm] * sin(θ) + x0
            smd.y[idx] = p.x[mm] * sin(θ) + p.y[mm] * cos(θ) + y0
        end
    end

    return smd
end

