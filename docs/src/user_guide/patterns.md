```@meta
CurrentModule = SMLMSim
```

# Pattern Types

SMLMSim includes several built-in pattern types for positioning fluorophores in both 2D and 3D space. This page explains how to use the different pattern types and how to create distributions of patterns.

## Pattern Hierarchy

All patterns in SMLMSim inherit from these abstract types:

- `Pattern`: Base type for all molecular patterns
  - `Pattern2D`: Base type for 2D patterns (x,y coordinates)
  - `Pattern3D`: Base type for 3D patterns (x,y,z coordinates)

## 2D Patterns

### Nmer2D

The `Nmer2D` pattern creates N molecules arranged symmetrically in a circle of specified diameter.

```julia
# Create an 8-molecule pattern with 100nm diameter (default)
nmer = Nmer2D()

# Create a custom pattern with 6 molecules and 200nm diameter
hexamer = Nmer2D(n=6, d=0.2)  # d is in microns
```

This is useful for simulating protein complexes, clusters, or other symmetric structures.

### Line2D

The `Line2D` pattern creates molecules arranged along a line segment with random positions.

```julia
# Create a line with default parameters (10 molecules/μm between (-1,0) and (1,0))
line = Line2D()

# Create a custom line (5 molecules/μm between (-2,0) and (2,0))
custom_line = Line2D(λ=5.0, endpoints=[(-2.0, 0.0), (2.0, 0.0)])
```

The density parameter λ (lambda) specifies the average number of molecules per micron along the line. The actual number of molecules is drawn from a Poisson distribution.

## 3D Patterns

### Nmer3D

The `Nmer3D` pattern is the 3D equivalent of `Nmer2D`, creating molecules in a circle at z=0.

```julia
# Create an 8-molecule pattern with 100nm diameter (default)
nmer3d = Nmer3D()

# Create a custom pattern with 6 molecules and 200nm diameter
hexamer3d = Nmer3D(n=6, d=0.2)
```

### Line3D

The `Line3D` pattern creates a 3D line with specified endpoints and density.

```julia
# Default 3D line along x-axis
line3d = Line3D()

# Custom 3D line with specified endpoints and density
custom_line3d = Line3D(
    λ=5.0,  # molecules per micron
    endpoints=[(-1.0, 0.0, -0.5), (1.0, 0.0, 0.5)]
)
```

## Creating Pattern Distributions

Individual patterns define the arrangement of a single group of molecules. To create realistic samples, you'll typically want to distribute many instances of a pattern throughout the field of view.

The `simulate()` function handles this automatically with the density parameter `ρ`:

```julia
# Simulate with 2 patterns per square micron
smld_true, smld_model, smld_noisy = simulate(
    ρ=2.0,
    pattern=Nmer2D(n=6, d=0.2)
)
```

Behind the scenes, this uses the `uniform2D()` or `uniform3D()` functions, which:

1. Generate a random number of patterns based on the density and field size
2. Place patterns randomly within the field of view
3. Apply random rotations to each pattern
4. Return the complete set of molecular coordinates

For advanced usage, you can access these functions directly:

```julia
# Create custom distribution of 2D patterns
field_x = 10.0  # μm
field_y = 10.0  # μm
pattern = Nmer2D(n=6, d=0.2)
density = 1.5   # patterns/μm²

x, y = uniform2D(density, pattern, field_x, field_y)

# Create custom distribution of 3D patterns
pattern3d = Nmer3D(n=6, d=0.2)
x, y, z = uniform3D(density, pattern3d, field_x, field_y, zrange=[-2.0, 2.0])
```

## Pattern Manipulation

### Rotation

You can rotate patterns before distributing them:

```julia
# Rotate a 2D pattern by 45 degrees
nmer = Nmer2D(n=8, d=0.1)
rotate!(nmer, π/4)

# Rotate a 3D pattern using Euler angles (ZYZ convention)
nmer3d = Nmer3D(n=8, d=0.1)
rotate!(nmer3d, π/4, π/6, π/3)  # α, β, γ angles in radians

# Rotate a 3D pattern using a rotation matrix
R = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]  # Z-axis rotation
rotate!(nmer3d, R)
```

### Custom Patterns

You can create your own pattern types by:

1. Define a new struct that inherits from `Pattern2D` or `Pattern3D`
2. Implement the necessary fields and constructors
3. Use your custom pattern with the existing simulation framework

Example:

```julia
# Create a custom 2D grid pattern
mutable struct Grid2D <: Pattern2D
    nx::Int  # number of columns
    ny::Int  # number of rows
    dx::Float64  # column spacing
    dy::Float64  # row spacing
    x::Vector{Float64}
    y::Vector{Float64}
end

function Grid2D(; nx=3, ny=3, dx=0.1, dy=0.1)
    n = nx * ny
    x = zeros(n)
    y = zeros(n)
    
    idx = 1
    for i in 1:nx, j in 1:ny
        x[idx] = (i - (nx+1)/2) * dx
        y[idx] = (j - (ny+1)/2) * dy
        idx += 1
    end
    
    return Grid2D(nx, ny, dx, dy, x, y)
end

# Use your custom pattern
grid = Grid2D(nx=4, ny=3, dx=0.1, dy=0.15)
smld_true, smld_model, smld_noisy = simulate(pattern=grid)
```

## Considerations for Realistic Simulations

When designing patterns for SMLM simulations, consider:

- **Labeling density**: Biological structures typically have incomplete labeling
- **Label size**: Fluorescent labels add ~5-10nm to the true structure size
- **Orientation**: Most structures have random orientations in samples
- **Size distribution**: Real structures often have size variations

Adjusting these parameters can help create more realistic simulations that better match experimental data.