# Pattern type definitions
abstract type Pattern end
abstract type Pattern2D <: Pattern end
abstract type Pattern3D <: Pattern end

#==========================================================================
2D Pattern Types
==========================================================================#

"""
    Nmer2D <: Pattern2D

N molecules symmetrically organized around a circle with diameter d    

# Fields
- `n`: Number of Points
- `d`: Diameter
- `x`: X positions
- `y`: Y positions
"""
mutable struct Nmer2D <: Pattern2D
    n::Int
    d::Float64
    x::Vector{Float64}
    y::Vector{Float64}
end

function Nmer2D(; n::Int=8, d::Float64=0.1)
    nmer = Nmer2D(n, d, zeros(n), zeros(n))
    for nn = 1:n
        θ = 2 * pi / n * (nn - 1)
        nmer.x[nn] = d / 2 * cos(θ)
        nmer.y[nn] = d / 2 * sin(θ)
    end
    return nmer
end

"""
    Nmer3D <: Pattern3D

N molecules symmetrically organized around a circle with diameter d at z=0    

# Fields
- `n`: Number of Points
- `d`: Diameter
- `x`: X positions
- `y`: Y positions
- `z`: Z positions
"""
mutable struct Nmer3D <: Pattern3D
    n::Int
    d::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end

function Nmer3D(; n::Int=8, d::Float64=0.1)
    nmer = Nmer3D(n, d, zeros(n), zeros(n), zeros(n))
    for nn = 1:n
        θ = 2 * pi / n * (nn - 1)
        nmer.x[nn] = d / 2 * cos(θ)
        nmer.y[nn] = d / 2 * sin(θ)
        nmer.z[nn] = 0.0
    end
    return nmer
end

"""
    Line2D <: Pattern2D

Points with uniform random distribution between 2 endpoints.    

# Fields
- `λ`: Linear molecule density
- `endpoints`: Vector of endpoint coordinates
- `n`: Number of Points
- `x`: X positions
- `y`: Y positions
"""
mutable struct Line2D <: Pattern2D
    n::Int
    x::Vector{Float64}
    y::Vector{Float64}
    λ::Float64
    endpoints::Vector{Tuple{Float64,Float64}}
end

function Line2D(; λ::Float64=10.0, endpoints=[(-1.0, 0.0), (1.0, 0.0)])
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
    Line3D <: Pattern3D

Points with uniform random distribution between 2 3D endpoints.    

# Fields
- `λ`: Linear molecule density
- `endpoints`: Vector of 3D endpoint coordinates
- `n`: Number of Points
- `x`: X positions
- `y`: Y positions
- `z`: Z positions
"""
mutable struct Line3D <: Pattern3D
    n::Int
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    λ::Float64
    endpoints::Vector{Tuple{Float64,Float64,Float64}}
end

function Line3D(; λ::Float64=10.0, endpoints=[(-1.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
    lx = (endpoints[2][1] - endpoints[1][1])
    ly = (endpoints[2][2] - endpoints[1][2])
    lz = (endpoints[2][3] - endpoints[1][3])
    l = sqrt(lx^2 + ly^2 + lz^2)

    pois = Poisson(λ * l)
    n = rand(pois)

    line = Line3D(n, zeros(n), zeros(n), zeros(n), λ, endpoints)
    for nn = 1:n
        d = l * rand()
        line.x[nn] = endpoints[1][1] + d / l * lx
        line.y[nn] = endpoints[1][2] + d / l * ly
        line.z[nn] = endpoints[1][3] + d / l * lz
    end
    return line
end

#==========================================================================
Pattern Generation Functions
==========================================================================#

"""
    uniform2D(ρ, p::Pattern2D, field_x::Float64, field_y::Float64)

Create coordinate arrays for randomly placed and rotated 2D patterns.

# Arguments
- `ρ`: Pattern density (patterns per square micron)
- `p`: Pattern to replicate
- `field_x`: Field width in microns
- `field_y`: Field height in microns

# Returns
Tuple{Vector{Float64}, Vector{Float64}}: (x, y) coordinates in microns
"""
function uniform2D(ρ::Float64, p::Pattern2D, field_x::Float64, field_y::Float64)
    # Generate random number of patterns
    npatterns = rand(Poisson(field_x * field_y * ρ))
    ntotal = npatterns * p.n

    # Initialize coordinate arrays
    x = Vector{Float64}(undef, ntotal)
    y = Vector{Float64}(undef, ntotal)
    
    idx = 1
    for nn = 1:npatterns
        θ = 2 * pi * rand()
        x0 = rand() * field_x
        y0 = rand() * field_y

        for mm = 1:p.n
            # Rotate and translate pattern points
            x[idx] = p.x[mm] * cos(θ) - p.y[mm] * sin(θ) + x0
            y[idx] = p.x[mm] * sin(θ) + p.y[mm] * cos(θ) + y0
            idx += 1
        end
    end

    return x, y
end

"""
    uniform3D(ρ, p::Pattern3D, field_x::Float64, field_y::Float64; 
             zrange::Vector{Float64}=[-1.0, 1.0])

Create coordinate arrays for randomly placed and rotated 3D patterns.

# Arguments
- `ρ`: Pattern density (patterns per square micron)
- `p`: Pattern to replicate
- `field_x`: Field width in microns
- `field_y`: Field height in microns
- `zrange`: [min_z, max_z] range in microns

# Returns
Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}: (x, y, z) coordinates
"""
function uniform3D(ρ::Float64, p::Pattern3D, field_x::Float64, field_y::Float64;
                  zrange::Vector{Float64}=Float64[-1.0, 1.0])
    # Generate random number of patterns (note: ρ is 2D density)
    npatterns = rand(Poisson(field_x * field_y * ρ))
    ntotal = npatterns * p.n

    # Initialize coordinate arrays
    x = Vector{Float64}(undef, ntotal)
    y = Vector{Float64}(undef, ntotal)
    z = Vector{Float64}(undef, ntotal)
    
    idx = 1
    for nn = 1:npatterns
        # Random position
        x0 = rand() * field_x
        y0 = rand() * field_y
        z0 = rand() * (zrange[2] - zrange[1]) + zrange[1]

        # Generate random 3D rotation using quaternions
        # This gives uniform rotation in 3D space
        u1, u2, u3 = rand(3)
        
        q0 = sqrt(1 - u1) * sin(2π * u2)
        q1 = sqrt(1 - u1) * cos(2π * u2)
        q2 = sqrt(u1) * sin(2π * u3)
        q3 = sqrt(u1) * cos(2π * u3)
        
        # Convert quaternion to rotation matrix
        R = [
            1-2*(q2^2 + q3^2)  2*(q1*q2 - q0*q3)  2*(q1*q3 + q0*q2);
            2*(q1*q2 + q0*q3)  1-2*(q1^2 + q3^2)  2*(q2*q3 - q0*q1);
            2*(q1*q3 - q0*q2)  2*(q2*q3 + q0*q1)  1-2*(q1^2 + q2^2)
        ]

        for mm = 1:p.n
            # Rotate and translate pattern points
            pos = R * [p.x[mm]; p.y[mm]; p.z[mm]]
            x[idx] = pos[1] + x0
            y[idx] = pos[2] + y0
            z[idx] = pos[3] + z0
            idx += 1
        end
    end

    return x, y, z
end

#==========================================================================
Pattern Rotation Functions
==========================================================================#

"""
    rotate!(p::Pattern2D, θ::Float64)

Rotate a 2D pattern by angle θ (in radians).
"""
function rotate!(p::Pattern2D, θ::Float64)
    for n in 1:p.n
        x_new = p.x[n] * cos(θ) - p.y[n] * sin(θ)
        y_new = p.x[n] * sin(θ) + p.y[n] * cos(θ)
        p.x[n] = x_new
        p.y[n] = y_new
    end
    return nothing
end

"""
    rotate!(p::Pattern3D, R::Matrix{Float64})

Rotate a 3D pattern by rotation matrix R.
"""
function rotate!(p::Pattern3D, R::Matrix{Float64})
    for n in 1:p.n
        pos = R * [p.x[n]; p.y[n]; p.z[n]]
        p.x[n] = pos[1]
        p.y[n] = pos[2]
        p.z[n] = pos[3]
    end
    return nothing
end

"""
    rotate!(p::Pattern3D, α::Float64, β::Float64, γ::Float64)

Rotate a 3D pattern by Euler angles α, β, γ (in radians).
Uses ZYZ convention.
"""
function rotate!(p::Pattern3D, α::Float64, β::Float64, γ::Float64)
    # ZYZ Euler angle rotation matrix
    R = [
        cos(α)*cos(β)*cos(γ)-sin(α)*sin(γ) -cos(α)*cos(β)*sin(γ)-sin(α)*cos(γ) cos(α)*sin(β);
        sin(α)*cos(β)*cos(γ)+cos(α)*sin(γ) -sin(α)*cos(β)*sin(γ)+cos(α)*cos(γ) sin(α)*sin(β);
        -sin(β)*cos(γ) sin(β)*sin(γ) cos(β)
    ]
    rotate!(p, R)
end