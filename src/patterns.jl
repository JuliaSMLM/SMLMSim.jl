#define patterned stuctures

"""
    Pattern

Abstract type for structured patterns of molecules    
"""
abstract type Pattern end


"""
    Nmer2D <: Pattern

N molecules symmetricaly organized around a circle with diameter d    

Nmer2D(;n::Int=8, d::AbstractFloat=.1)

# Fields
- 'n': Numbor of Points = 1
- 'd': Diameter
- 'x': X position
- 'y': Y position

"""
mutable struct Nmer2D <: Pattern
    n::Int
    d::AbstractFloat
    x::Vector{AbstractFloat}
    y::Vector{AbstractFloat}
end
function Nmer2D(;n::Int=8, d::AbstractFloat=.1)

    nmer = Nmer2D(n, d, zeros(n), zeros(n))
    for nn = 1:n
        θ = 2 * pi / n * (nn - 1)
        nmer.x[nn] = d / 2 * cos(θ)
        nmer.y[nn] = d / 2 * sin(θ)
    end
    return nmer
end



"""
    Point2D <: Pattern

    A single 2D point.    

Point2D() = new(1, [0.0], [0.0])

# Fields
- 'n': Numbor of Points = 1
- 'x': X position
- 'y': Y position

"""
mutable struct Point2D <: Pattern
    n::Int
    x::Vector{AbstractFloat}
    y::Vector{AbstractFloat}
    Point2D() = new(1, [0.0], [0.0])
end


"""
    Line2D <: Pattern

Points with uniform random distribution between 2 endpoints.    

Line2D(;λ::AbstractFloat=10.0, endpoints=[(-1.0,0.0),(1.0,0.0)])

# Fields
- `λ`: linear molecule density
- 'endpoints': Vector of Tuple 
- 'n': Numbor of Points = 1
- 'x': X position
- 'y': Y position

"""
mutable struct Line2D <: Pattern
    n::Int
    x::Vector{AbstractFloat}
    y::Vector{AbstractFloat}
    λ::AbstractFloat
    endpoints::Vector{Tuple{<:AbstractFloat,<:AbstractFloat}}
end
function Line2D(;λ::AbstractFloat=10.0, endpoints=[(-1.0,0.0),(1.0,0.0)])
    
    lx=(endpoints[2][1]-endpoints[1][1])
    ly=(endpoints[2][2]-endpoints[2][1])
    l=sqrt( lx^2+ ly^2)

    pois=Poisson(λ*l)
    n=rand(pois)

    line = Line2D(n, zeros(n), zeros(n),λ,endpoints)
    for nn = 1:n
        d=l*rand()
        line.x[nn] = endpoints[1][1] + d/l*lx
        line.y[nn] = endpoints[1][2] + d/l*ly
    end
    return line
end





"""
    function uniformPattern2D(ρ,xsize::AbstractFloat,ysize::AbstractFloat, p::Pattern)

Create positions of molecules from uniformly randomly placed and rotated patterns
"""
function uniform2D(ρ, p::Pattern, xsize::Real, ysize::Real)

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

## Pattern Rotation

"""
    rotate!(p::Pattern,θ::AbstractFloat)

Rotate a Pattern in 2D by \\theta.


Both molecule positions and reference positions are rotated (e.g. endpoints of a line)
"""
function rotate!(p::Pattern,θ::AbstractFloat)
end

"""
Rotate a Pattern in 3D by the improper Euler angles [\\alpha \\beta \\gamma].
"""
function rotate!(p::Pattern,α::Real, β::Real, γ::Real)
end


"""
Rotate a Pattern in 3D by premultiplying with the rotation matrix `r`
"""
function rotate!(p::Pattern,r::Array{AbstractFloat})
end

function rotate(x::Real,y::Real,θ::Real)
    return (x*cos(θ)-y*sin(θ), x*sin(θ)+y*cos(θ))
end

function rotate(x::Real,y::Real,z::Real,r::Matrix{<:Real})
    out = r*[x y z]'
    return (out[1],out[2],out[3])       
end

function rotate(x::Real,y::Real,z::Real,α::Real, β::Real, γ::Real)
    r = [
        cos(β)*cos(γ) sin(α)*sin(β)*cos(γ)-cos(α)*sin(γ) cos(α)*sin(β)*cos(γ)+sin(α)*sin(γ)
        cos(β)*sin(γ) sin(α)*sin(β)*sin(γ)+cos(α)*cos(γ) cos(α)*sin(β)*sin(γ)-sin(α)*cos(γ)
        -sin(β) sin(α)*cos(β) cos(α)*cos(β)
        ]
    return rotate(x,y,z,r)
end

function rotate!(p::Point2D)
    

end

