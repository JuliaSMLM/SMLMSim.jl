```@meta
CurrentModule = SMLMSim
```

# Performance Tips

This page provides strategies for optimizing SMLMSim simulations, particularly for large-scale or computationally intensive scenarios.

## Simulation Scale Considerations

SMLM simulations can become computationally intensive as you increase:
- Number of emitters
- Number of frames
- Number of datasets
- Complexity of patterns or dynamics

Here are guidelines for managing simulation scale:

| Parameter | Small | Medium | Large |
|-----------|-------|--------|-------|
| Emitter density (ρ) | 0.1-1.0/μm² | 1.0-10.0/μm² | >10.0/μm² |
| Field size | <50×50 μm | 50-500×50-500 μm | >500×500 μm |
| Frames | <1,000 | 1,000-10,000 | >10,000 |
| Datasets | 1-10 | 10-100 | >100 |

## Memory Management

### Controlling Memory Usage

Large simulations can consume significant memory, especially with high emitter densities and frame counts:

```julia
# Memory-intensive simulation (use caution)
smld_true, smld_model, smld_noisy = simulate(
    ρ=5.0,                  # 5 patterns per μm²
    nframes=10000,          # 10,000 frames
    ndatasets=10,           # 10 independent datasets
    camera=IdealCamera(1:1024, 1:1024, 0.1)  # Large field of view
)
```

To reduce memory usage:

1. **Decrease field of view** when possible
   ```julia
   camera = IdealCamera(128, 128, 0.1)  # Smaller field = fewer emitters
   ```

2. **Process data incrementally** by simulating one dataset at a time
   ```julia
   results = []
   for i in 1:10
       _, _, smld_noisy = simulate(ndatasets=1, ...)
       # Process dataset immediately
       processed_data = analyze_dataset(smld_noisy)
       push!(results, processed_data)
       # Allow garbage collection between iterations
       GC.gc()
   end
   ```

3. **Optimize photophysics parameters** to reduce the number of active emitters per frame
   ```julia
   function run_sim(; use_generic_fluor=true)
       if use_generic_fluor
           # Use positional constructor
           molecule = GenericFluor(1e4, [-50.0 50.0; 0.5 -0.5]) # γ=1e4, k_off=50, k_on=0.5
       else
           # Use 2-state keyword constructor
           molecule = GenericFluor(; photons=1e4, k_off=50.0, k_on=0.5)
       end
   end
   ```

4. **Avoid storing intermediate results** when possible
   ```julia
   # If you only need the final noisy result
   _, _, smld_noisy = simulate(...)
   ```

## Computational Performance

### Parallelization

SMLMSim uses Julia's multi-threading capabilities for some computations. Enable multiple threads when running Julia:

```bash
julia --threads=auto  # Use all available cores
```

You can check the number of threads available:

```julia
Threads.nthreads()
```

### Precompilation

To reduce runtime, precompile the package and key functions:

```julia
# Warmup run with minimal parameters
camera = IdealCamera(1:32, 1:32, 0.1)
_ = simulate(ρ=0.1, nframes=10, camera=camera)

# Now run your actual simulation
smld_true, smld_model, smld_noisy = simulate(...)
```

### GPU Acceleration

While SMLMSim doesn't currently offer direct GPU support, you can use GPU-accelerated libraries for post-processing:

```julia
using CUDA  # Requires CUDA.jl
using CuArrays

# Convert a large dataset to GPU arrays for faster processing
x_gpu = CuArray([e.x for e in smld_noisy.emitters])
y_gpu = CuArray([e.y for e in smld_noisy.emitters])

# Process on GPU
# ...

# Return to CPU when finished
x_cpu = Array(x_gpu)
y_cpu = Array(y_gpu)
```

## Optimizing Specific Simulations

### Interaction-Diffusion Simulations

The Smoluchowski diffusion simulations can be particularly computationally intensive:

1. **Optimize time step**: Larger time steps (`dt`) reduce computation time but decrease accuracy
   ```julia
   params = SmoluchowskiParams(
       density = 0.5,
       dt = 0.05,           # Larger time step = faster but less accurate
       t_max = 10.0
   )
   ```

2. **Limit box size and molecule count** for fast prototyping
   ```julia
   params = SmoluchowskiParams(
       density = 0.2,       # Lower density = fewer molecules
       box_size = 5.0,      # Smaller box = fewer molecules
       t_max = 5.0          # Shorter simulation time
   )
   ```

3. **Optimize integration steps** when generating microscope images
   ```julia
   # Fewer integration steps = faster rendering
   images = gen_image_sequence(psf, systems, frame_integration=5)
   ```

### Large 3D Simulations

For 3D simulations with many emitters:

1. **Optimize axial range** to include only the region of interest
   ```julia
   # More focused z-range
   smld_true, smld_model, smld_noisy = simulate(
       pattern=Nmer3D(),
       zrange=[-0.5, 0.5]  # More limited axial range
   )
   ```

2. **Use appropriate density** for the study
   ```julia
   # Lower density for 3D simulations
   smld_true, smld_model, smld_noisy = simulate(
       pattern=Nmer3D(),
       ρ=0.3  # Patterns per μm²
   )
   ```

## Benchmarking Your Simulations

You can benchmark performance to identify bottlenecks:

```julia
using BenchmarkTools

# Benchmark pattern generation
@btime uniform2D(1.0, Nmer2D(), 10.0, 10.0);

# Benchmark kinetic model
@btime kinetic_model($smld_true, $fluor, 1000, 50.0);

# Benchmark noise addition
@btime noise($smld_model, 0.13);
```

## Julia-Specific Optimizations

To get the best performance from SMLMSim:

1. **Avoid global variables** in performance-critical code
   ```julia
   # Good: Pass parameters directly
   function analyze_simulation(smld_noisy, threshold)
       # ...
   end
   
   # Bad: Reference global variables
   threshold = 100
   function analyze_simulation(smld_noisy)
       # Using global threshold
       # ...
   end
   ```

2. **Use type stability** in custom analysis functions
   ```julia
   # Good: Type-stable function
   function calculate_distances(emitters::Vector{<:Emitter2D{Float64}})
       # ...
   end
   
   # Bad: Type-unstable function
   function calculate_distances(emitters)
       # ...
   end
   ```

3. **Preallocate arrays** for accumulating results
   ```julia
   # Good: Preallocate
   function trajectory_analysis(emitters, n_frames)
       positions = Vector{Tuple{Float64, Float64}}(undef, n_frames)
       # ...
       return positions
   end
   
   # Bad: Grow array dynamically
   function trajectory_analysis(emitters, n_frames)
       positions = []
       for i in 1:n_frames
           # push! to positions
       end
       return positions
   end
   ```

By following these guidelines, you can significantly improve the performance of your SMLM simulations, especially for large-scale or complex scenarios.