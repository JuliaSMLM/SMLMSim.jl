# SMLMSim.jl Development Guide

## Build, Test & Run Commands
- Install dependencies: `julia --project -e "using Pkg; Pkg.instantiate()"`
- Run all tests: `julia --project -e "using Pkg; Pkg.test()"`
- Run specific test: `julia --project -e "using Pkg; Pkg.test(\"SMLMSim\", test_args=[\"Patterns\"])"`
- Build docs: `julia --project=docs/ docs/make.jl`
- Run benchmark: `julia --project=dev/ dev/benchmark.jl`

## Code Style Guidelines
- Use physical units consistently (μm for space, seconds for time)
- Follow Julia naming conventions: lowercase for functions, CamelCase for types
- Organize imports: standard library first, then external packages
- Place exports at module level, grouped by functionality
- Type parameters should be explicit for functions with multiple dispatch
- Error handling: use descriptive errors with parameter names and constraints
- Documentation: include docstrings with examples and physical units
- Re-export SMLMData types rather than redefining them
- Follow existing patterns for defining new molecule or pattern types
- Tests should use the `@test` macro with appropriate tolerances for floating-point comparisons

## Module Structure
Functions should be organized in appropriate submodules: core, static, or diffusion.