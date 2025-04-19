# Pattern type definitions
"""
    Pattern

Abstract type for all molecular spatial patterns.
"""
abstract type Pattern end

"""
    Pattern2D <: Pattern

Abstract type for 2D molecular spatial patterns.
"""
abstract type Pattern2D <: Pattern end

"""
    Pattern3D <: Pattern

Abstract type for 3D molecular spatial patterns.
"""
abstract type Pattern3D <: Pattern end

#==========================================================================
2D Pattern Types
==========================================================================#

"""
    Nmer2D <: Pattern2D

N molecules symmetrically organized around a circle with diameter d.

# Fields
- `n::Int`: Number of molecules in the pattern
- `d::Float64`: Diameter of the circle in microns
- `x::Vector{Float64}`: X positions of molecules in microns
- `y::Vector{Float64}`: Y positions of molecules in microns

# Examples
```julia
# Create an 8-molecule pattern with 100nm diameter
nmer = Nmer2D()

# Create a custom pattern with 6 molecules and 200nm diameter
nmer = Nmer2D(; n=6, d=0.2)
```
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

function Base.show(io::IO, nmer::Nmer2D)
    print(io, "Nmer2D(n=$(nmer.n), d=$(nmer.d) μm)")
end

function Base.show(io::IO, ::MIME"text/plain", nmer::Nmer2D)
    println(io, "Nmer2D pattern:")
    println(io, "  Number of molecules (n) = $(nmer.n)")
    println(io, "  Diameter (d) = $(nmer.d) μm ($(nmer.d*1000) nm)")
    
    if nmer.n <= 10  # Only show coordinates for small patterns
        println(io, "  Coordinates (μm):")
        for i in 1:nmer.n
            print(io, "    [$(i)]: (")
            print(io, @sprintf("%.3f", nmer.x[i]))
            print(io, ", ")
            print(io, @sprintf("%.3f", nmer.y[i]))
            println(io, ")")
        end
    end
end

"""
    Nmer3D <: Pattern3D

N molecules symmetrically organized around a circle with diameter d at z=0.

# Fields
- `n::Int`: Number of molecules in the pattern
- `d::Float64`: Diameter of the circle in microns
- `x::Vector{Float64}`: X positions of molecules in microns
- `y::Vector{Float64}`: Y positions of molecules in microns
- `z::Vector{Float64}`: Z positions of molecules in microns

# Examples
```julia
# Create an 8-molecule pattern with 100nm diameter
nmer = Nmer3D()

# Create a custom pattern with 6 molecules and 200nm diameter
nmer = Nmer3D(; n=6, d=0.2)
```
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

function Base.show(io::IO, nmer::Nmer3D)
    print(io, "Nmer3D(n=$(nmer.n), d=$(nmer.d) μm)")
end

function Base.show(io::IO, ::MIME"text/plain", nmer::Nmer3D)
    println(io, "Nmer3D pattern:")
    println(io, "  Number of molecules (n) = $(nmer.n)")
    println(io, "  Diameter (d) = $(nmer.d) μm ($(nmer.d*1000) nm)")
    
    if nmer.n <= 10  # Only show coordinates for small patterns
        println(io, "  Coordinates (μm):")
        for i in 1:nmer.n
            print(io, "    [$(i)]: (")
            print(io, @sprintf("%.3f", nmer.x[i]))
            print(io, ", ")
            print(io, @sprintf("%.3f", nmer.y[i]))
            print(io, ", ")
            print(io, @sprintf("%.3f", nmer.z[i]))
            println(io, ")")
        end
    end
end

"""
    Line2D <: Pattern2D

Points with uniform random distribution between two endpoints.

# Fields
- `λ::Float64`: Linear molecule density (molecules per micron)
- `endpoints::Vector{Tuple{Float64,Float64}}`: Vector of endpoint coordinates
- `n::Int`: Number of molecules in the pattern
- `x::Vector{Float64}`: X positions of molecules in microns
- `y::Vector{Float64}`: Y positions of molecules in microns

# Examples
```julia
# Create a line with default parameters
line = Line2D()

# Create a custom line
line = Line2D(; λ=5.0, endpoints=[(-2.0, 0.0), (2.0, 0.0)])
```
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

function Base.show(io::IO, line::Line2D)
    p1 = line.endpoints[1]
    p2 = line.endpoints[2]
    length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
    print(io, "Line2D(n=$(line.n), λ=$(line.λ)/μm, length=$(round(length, digits=2)) μm)")
end

function Base.show(io::IO, ::MIME"text/plain", line::Line2D)
    p1 = line.endpoints[1]
    p2 = line.endpoints[2]
    length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2)
    
    println(io, "Line2D pattern:")
    println(io, "  Number of molecules (n) = $(line.n)")
    println(io, "  Linear density (λ) = $(line.λ) molecules/μm")
    println(io, "  Length = $(round(length, digits=3)) μm")
    println(io, "  Endpoints:")
    println(io, "    Start: ($(p1[1]), $(p1[2])) μm")
    println(io, "    End:   ($(p2[1]), $(p2[2])) μm")
    
    if line.n <= 10  # Only show coordinates for small patterns
        println(io, "  Molecule coordinates (μm):")
        for i in 1:line.n
            print(io, "    [$(i)]: (")
            print(io, @sprintf("%.3f", line.x[i]))
            print(io, ", ")
            print(io, @sprintf("%.3f", line.y[i]))
            println(io, ")")
        end
    end
end

"""
    Line3D <: Pattern3D

Points with uniform random distribution between two 3D endpoints.

# Fields
- `λ::Float64`: Linear molecule density (molecules per micron)
- `endpoints::Vector{Tuple{Float64,Float64,Float64}}`: Vector of 3D endpoint coordinates
- `n::Int`: Number of molecules in the pattern
- `x::Vector{Float64}`: X positions of molecules in microns
- `y::Vector{Float64}`: Y positions of molecules in microns
- `z::Vector{Float64}`: Z positions of molecules in microns

# Examples
```julia
# Create a line with default parameters
line = Line3D()

# Create a custom 3D line
line = Line3D(; λ=5.0, endpoints=[(-1.0, 0.0, -0.5), (1.0, 0.0, 0.5)])
```
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

function Base.show(io::IO, line::Line3D)
    p1 = line.endpoints[1]
    p2 = line.endpoints[2]
    length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2 + (p2[3] - p1[3])^2)
    print(io, "Line3D(n=$(line.n), λ=$(line.λ)/μm, length=$(round(length, digits=2)) μm)")
end

function Base.show(io::IO, ::MIME"text/plain", line::Line3D)
    p1 = line.endpoints[1]
    p2 = line.endpoints[2]
    length = sqrt((p2[1] - p1[1])^2 + (p2[2] - p1[2])^2 + (p2[3] - p1[3])^2)
    
    println(io, "Line3D pattern:")
    println(io, "  Number of molecules (n) = $(line.n)")
    println(io, "  Linear density (λ) = $(line.λ) molecules/μm")
    println(io, "  Length = $(round(length, digits=3)) μm")
    println(io, "  Endpoints:")
    println(io, "    Start: ($(p1[1]), $(p1[2]), $(p1[3])) μm")
    println(io, "    End:   ($(p2[1]), $(p2[2]), $(p2[3])) μm")
    
    if line.n <= 10  # Only show coordinates for small patterns
        println(io, "  Molecule coordinates (μm):")
        for i in 1:line.n
            print(io, "    [$(i)]: (")
            print(io, @sprintf("%.3f", line.x[i]))
            print(io, ", ")
            print(io, @sprintf("%.3f", line.y[i]))
            print(io, ", ")
            print(io, @sprintf("%.3f", line.z[i]))
            println(io, ")")
        end
    end
end

#==========================================================================
Pattern Generation Functions
==========================================================================#

"""
    uniform2D(ρ::Float64, p::Pattern2D, field_x::Float64, field_y::Float64)

Create coordinate arrays for randomly placed and rotated 2D patterns.

# Arguments
- `ρ::Float64`: Pattern density (patterns per square micron)
- `p::Pattern2D`: Pattern to replicate
- `field_x::Float64`: Field width in microns
- `field_y::Float64`: Field height in microns

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: (x, y) coordinates in microns

# Example
```julia
# Generate coordinates for randomly placed Nmer2D patterns
nmer = Nmer2D(; n=6, d=0.2)
x, y = uniform2D(1.0, nmer, 10.0, 10.0)
```
"""
function uniform2D(ρ::T, p::Pattern2D, field_x::T, field_y::T) where T <: AbstractFloat
    # Generate random number of patterns
    npatterns = rand(Poisson(field_x * field_y * ρ))
    ntotal = npatterns * p.n

    # Initialize coordinate arrays
    x = Vector{T}(undef, ntotal)
    y = Vector{T}(undef, ntotal)
    
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
    uniform3D(ρ::Float64, p::Pattern3D, field_x::Float64, field_y::Float64; 
             zrange::Vector{Float64}=[-1.0, 1.0])

Create coordinate arrays for randomly placed and rotated 3D patterns.

# Arguments
- `ρ::Float64`: Pattern density (patterns per square micron)
- `p::Pattern3D`: Pattern to replicate
- `field_x::Float64`: Field width in microns
- `field_y::Float64`: Field height in microns
- `zrange::Vector{Float64}=[-1.0, 1.0]`: [min_z, max_z] range in microns

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}`: (x, y, z) coordinates

# Example
```julia
# Generate coordinates for randomly placed Nmer3D patterns
nmer = Nmer3D(; n=6, d=0.2)
x, y, z = uniform3D(1.0, nmer, 10.0, 10.0; zrange=[-2.0, 2.0])
```
"""
function uniform3D(ρ::T, p::Pattern3D, field_x::T, field_y::T;
                  zrange::Vector{T}=T[-1.0, 1.0]) where T <: AbstractFloat
    # Input validation
    if length(zrange) != 2 || zrange[1] >= zrange[2]
        throw(ArgumentError("zrange must be a vector of two values [min_z, max_z] where min_z < max_z"))
    end
    
    # Generate random number of patterns (note: ρ is 2D density)
    npatterns = rand(Poisson(field_x * field_y * ρ))
    ntotal = npatterns * p.n

    # Initialize coordinate arrays
    x = Vector{T}(undef, ntotal)
    y = Vector{T}(undef, ntotal)
    z = Vector{T}(undef, ntotal)
    
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

# Arguments
- `p::Pattern2D`: Pattern to rotate
- `θ::Float64`: Rotation angle in radians

# Example
```julia
nmer = Nmer2D()
rotate!(nmer, π/4)  # Rotate 45 degrees
```
"""
function rotate!(p::Pattern2D, θ::T) where T <: AbstractFloat
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

# Arguments
- `p::Pattern3D`: Pattern to rotate
- `R::Matrix{Float64}`: 3x3 rotation matrix

# Example
```julia
nmer = Nmer3D()
# Create a rotation matrix for 90 degrees around z-axis
θ = π/2
R = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
rotate!(nmer, R)
```
"""
function rotate!(p::Pattern3D, R::Matrix{T}) where T <: AbstractFloat
    if size(R) != (3, 3)
        throw(ArgumentError("Rotation matrix must be 3×3"))
    end
    
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

# Arguments
- `p::Pattern3D`: Pattern to rotate
- `α::Float64`: First rotation angle (around Z axis)
- `β::Float64`: Second rotation angle (around Y' axis)
- `γ::Float64`: Third rotation angle (around Z'' axis)

# Example
```julia
nmer = Nmer3D()
rotate!(nmer, π/4, π/6, π/3)
```
"""
function rotate!(p::Pattern3D, α::T, β::T, γ::T) where T <: AbstractFloat
    # ZYZ Euler angle rotation matrix
    R = [
        cos(α)*cos(β)*cos(γ)-sin(α)*sin(γ) -cos(α)*cos(β)*sin(γ)-sin(α)*cos(γ) cos(α)*sin(β);
        sin(α)*cos(β)*cos(γ)+cos(α)*sin(γ) -sin(α)*cos(β)*sin(γ)+cos(α)*cos(γ) sin(α)*sin(β);
        -sin(β)*cos(γ) sin(β)*sin(γ) cos(β)
    ]
    rotate!(p, R)
end