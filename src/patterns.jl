#define patterned stuctures

"""
    Pattern

Abstract type for structured patterns of molecules    
"""
abstract type Pattern end

"""
    Pattern2D

Abstract type for structured patterns of molecules    
"""
abstract type Pattern2D <: Pattern end

"""
    Pattern3D

Abstract type for structured patterns of molecules    
"""
abstract type Pattern3D <: Pattern end




"""
    Nmer2D <: Pattern2D

N molecules symmetricaly organized around a circle with diameter d    

Nmer2D(;n::Int=8, d::AbstractFloat=.1)

# Fields
- 'n': Numbor of Points = 1
- 'd': Diameter
- 'x': X position
- 'y': Y position

"""
mutable struct Nmer2D <: Pattern2D
    n::Int
    d::AbstractFloat
    x::Vector{AbstractFloat}
    y::Vector{AbstractFloat}
end
function Nmer2D(;  n::Int=8, d::AbstractFloat=00.1)

    nmer = Nmer2D(n, d, zeros(n), zeros(n))
    for nn = 1:n
        θ = 2 * pi / n * (nn - 1)
        nmer.x[nn] = d / 2 * cos(θ)
        nmer.y[nn] = d / 2 * sin(θ)
    end
    return nmer
end



"""
    Point2D <: Pattern2D

    A single 2D point.    

Point2D() = new(1, [0.0], [0.0])

# Fields
- 'n': Numbor of Points = 1
- 'x': X position
- 'y': Y position

"""
mutable struct Point2D <: Pattern2D
    n::Int
    x::Vector{AbstractFloat}
    y::Vector{AbstractFloat}
    Point2D() = new(1, [0.0], [0.0])
end

"""
    Point3D <: Pattern3D

    A single 3D point.    

Point3D() = new(1, [0.0], [0.0],[0.0])

# Fields
- 'n': Numbor of Points = 1
- 'x': X position
- 'y': Y position
- 'z': Z position
"""
mutable struct Point3D <: Pattern3D
    n::Int
    x::Vector{AbstractFloat}
    y::Vector{AbstractFloat}
    z::Vector{AbstractFloat}
    Point3D() = new(1, [0.0], [0.0], [0.0])
end


"""
    Line2D <: Pattern2D

Points with uniform random distribution between 2 endpoints.    

Line2D(;λ::AbstractFloat=10.0, endpoints=[(-1.0,0.0),(1.0,0.0)])

# Fields
- `λ`: linear molecule density
- 'endpoints': Vector of Tuple 
- 'n': Numbor of Points = 1
- 'x': X position
- 'y': Y position

"""
mutable struct Line2D <: Pattern2D
    n::Int
    x::Vector{AbstractFloat}
    y::Vector{AbstractFloat}
    λ::AbstractFloat
    endpoints::Vector{Tuple{<:AbstractFloat,<:AbstractFloat}}
end
function Line2D(; λ::AbstractFloat=10.0, endpoints=[(-1.0, 0.0), (1.0, 0.0)])

    lx = (endpoints[2][1] - endpoints[1][1])
    ly = (endpoints[2][2] - endpoints[1][2])
    l = sqrt(lx^2 + ly^2)

    pois = Poisson(λ * l)
    n = rand(pois)

    line = Line2D(n, zeros(n), zeros(n), λ, endpoints)
    for nn = 1:n
        d = l * rand()
        line.x[nn] = endpoints[1][1] + d / l * lx
        line.y[nn] = endpoints[1][2] + d / l * ly
    end
    return line
end





"""
    uniform2D(ρ, p::Pattern2D, xsize::AbstractFloat,ysize::AbstractFloat)

Create positions of molecules from uniformly randomly placed and rotated patterns.
"""
function uniform2D(ρ, p::Pattern2D, xsize::Real, ysize::Real)

    npatterns = rand(Poisson(xsize * ysize * ρ))
    ntotal = npatterns * p.n

    #make smd 
    smd = SMLMData.SMLD2D(ntotal)
    smd.datasize = Int.(ceil.([ysize; xsize]))
    for nn = 1:npatterns
        θ = 2 * pi * rand()
        x0 = rand() * xsize
        y0 = rand() * ysize

        for mm = 1:p.n
            idx = (p.n) * (nn - 1) + mm
            smd.x[idx] = p.x[mm] * cos(θ) - p.y[mm] * sin(θ) + x0
            smd.y[idx] = p.x[mm] * sin(θ) + p.y[mm] * cos(θ) + y0
        end
    end

    return smd.y, smd.x
end




"""
    uniform2D(p::Vector{Pattern2D}, xsize::Real, ysize::Real)  

Randomly place and rotate the input patterns.
"""
function uniform2D(p::Vector{Pattern2D},  xsize::Real, ysize::Real)

    npatterns = length(p)

    ntotal = 0
    for n in 1:npatterns
        ntotal += p[n].n
    end

    #make smd 
    smd = SMLMData.SMLD2D(ntotal)
    smd.datasize = Int.(ceil.([ysize; xsize]))
    idx=1
    for nn = 1:npatterns
        θ = 2 * pi * rand()
        x0 = rand() * xsize
        y0 = rand() * ysize
        for mm = 1:p.n
            smd.x[idx] = p[nn].x[mm] * cos(θ) - p[nn].y[mm] * sin(θ) + x0
            smd.y[idx] = p[nn].x[mm] * sin(θ) + p[nn].y[mm] * cos(θ) + y0
            idx += 1
        end
    end

    return smd.y, smd.x
end

"""
    uniform2D(p::Vector{Pattern}, xsize::Real, ysize::Real,θ::Real)  

Randomly place the input patterns with fixed rotation \\theta.
"""
function uniform2D(p::Vector{<:Pattern2D},  xsize::Real, ysize::Real, θ::Real)

    npatterns = length(p)

    ntotal = 0
    for n in 1:npatterns
        ntotal += p[n].n
    end

    #make smd 
    smd = SMLMData.SMLD2D(ntotal)
    smd.datasize = Int.(ceil.([ysize; xsize]))
    idx=1
    for nn = 1:npatterns
        x0 = rand() * xsize
        y0 = rand() * ysize
        for mm = 1:p[nn].n
            smd.x[idx] = p[nn].x[mm] * cos(θ) - p[nn].y[mm] * sin(θ) + x0
            smd.y[idx] = p[nn].x[mm] * sin(θ) + p[nn].y[mm] * cos(θ) + y0
            idx += 1
        end
    end

    return smd.y, smd.x
end

"""
    place2D(p::Vector{Pattern}, xsize::Real, ysize::Real)  

Place the input patterns and return SMLD2D.
"""
function place2D(p::Vector{<:Pattern2D},  xsize::Real, ysize::Real)

    npatterns = length(p)

    ntotal = 0
    for n in 1:npatterns
        ntotal += p[n].n
    end

    #make smd 
    smd = SMLMData.SMLD2D(ntotal)
    smd.datasize = Int.(ceil.([ysize; xsize]))
    idx=1
    for nn in 1:npatterns
        for mm in 1:p[nn].n
            smd.x[idx] = p[nn].x[mm] 
            smd.y[idx] = p[nn].y[mm] 
            idx += 1
        end
    end

   
    return smd.y, smd.x
end





"""
    function uniform3D(ρ,p::Pattern, xsize::Real,ysize::Real; zrange::Vector{<:Real}=[-1.0,1.0])

Create positions of molecules from uniformly randomly placed and rotated patterns.

!!! Note
`ρ` is 2D density. 3D density is `ρ/(zrange[2]-zrange[1])`. 

"""
function uniform3D(ρ, p::Pattern3D, xsize::Real, ysize::Real; zrange::Vector{<:Real}=[-1.0, 1.0])

    npatterns = rand(Poisson(xsize * ysize * ρ))
    ntotal = npatterns * p.n

    #make smd 
    smd = SMLMData.SMLD3D(ntotal)
    smd.datasize = Int.(ceil.([ysize; xsize]))
    for nn = 1:npatterns

        x0 = rand() * xsize
        y0 = rand() * ysize
        z0 = rand() * (zrange[2] - zrange[1]) + zrange[1]

        #Transformation that gives uniform rotation in 3D
        # based on J. Avro 1992
        x1 = rand()
        x2 = rand()
        x3 = rand()

        r = [
            cos(2 * pi * x1) -sin(2 * pi * x1) 0
            -sin(2 * pi * x1) cos(2 * pi * x1) 0
            0 0 1
        ]

        v = [
            cos(2 * pi * x2) * sqrt(x3)
            sin(2 * pi * x2) * sqrt(x3)
            sqrt(1 - x3)
        ]

        h = Diagonal([1.0, 1.0, 1.0]) - 2 * v * v'

        m = -r * h

        for mm in 1:p.n
            idx = (p.n) * (nn - 1) + mm

            xyz = [p.x[mm], p.y[mm], p.z[mm]]
            xyz_prime = m * xyz

            smd.x[idx] = xyz_prime[1] + x0
            smd.y[idx] = xyz_prime[2] + y0
            smd.z[idx] = xyz_prime[3] + z0
        end
    end

    return smd.y, smd.x, smd.z
end

## Pattern Rotation

"""
    rotate!(p::Pattern,θ::Real)

Rotate a Pattern in 2D by \\theta radians.


Both molecule positions and reference positions are rotated (e.g. endpoints of a line)
"""
function rotate!(p::Pattern, θ::Real)
end

"""
Rotate a Pattern in 3D by the improper Euler angles [\\alpha \\beta \\gamma] (radians).
"""
function rotate!(p::Pattern, α::Real, β::Real, γ::Real)
end


"""
Rotate a Pattern in 3D by premultiplying with the rotation matrix `r`.
"""
function rotate!(p::Pattern, r::Array{Real})
end

function rotate(x::Real, y::Real, θ::Real)
    return (x * cos(θ) - y * sin(θ), x * sin(θ) + y * cos(θ))
end

function rotate(x::Real, y::Real, z::Real, r::Matrix{<:Real})
    out = r * [x y z]'
    return (out[1], out[2], out[3])
end

function rotate(x::Real, y::Real, z::Real, α::Real, β::Real, γ::Real)
    r = [
        cos(β)*cos(γ) sin(α)*sin(β)*cos(γ)-cos(α)*sin(γ) cos(α)*sin(β)*cos(γ)+sin(α)*sin(γ)
        cos(β)*sin(γ) sin(α)*sin(β)*sin(γ)+cos(α)*cos(γ) cos(α)*sin(β)*sin(γ)-sin(α)*cos(γ)
        -sin(β) sin(α)*cos(β) cos(α)*cos(β)
    ]
    return rotate(x, y, z, r)
end

function rotate!(p::Point2D, θ::Real)
    for n in 1:p.n
        (p.x[n], p.y[n]) = rotate(p.x[n], p.y[n], θ)
    end
    return nothing
end

function rotate!(p::Line2D, θ::Real)
    for n in 1:p.n
        (p.x[n], p.y[n]) = rotate(p.x[n], p.y[n], θ)
    end
    p.endpoints[1] = rotate(p.endpoints[1][1], p.endpoints[1][2], θ)
    p.endpoints[2] = rotate(p.endpoints[2][1], p.endpoints[2][2], θ)
    return nothing
end



