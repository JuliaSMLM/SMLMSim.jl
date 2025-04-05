using SMLMSim
using SMLMData
using MicroscopePSFs

println("Testing SMLMSim integration with SMLMData v0.2.2...")

# Test 1: Create a basic simulation
println("Test 1: Creating a basic simulation...")
camera = IdealCamera(1:128, 1:128, 0.1)  # 128×128 pixels, 100nm pixels
params = StaticSMLMParams(ρ=1.0)  # Default parameters
println("  Creating simulation parameters: ", typeof(params))

# Test 2: Test CTMC functionality
println("Test 2: Testing CTMC functionality...")
q = [0.0 1.0; 0.1 0.0]
ctmc = CTMC(q, 10.0, 1)
state = get_state(ctmc, 0.5)
next_state, next_time = get_next(ctmc, 0.5)
println("  CTMC state at t=0.5: ", state)
println("  Next state after t=0.5: ", next_state, " at time ", next_time)

# Test 3: Create a pattern
println("Test 3: Testing pattern creation...")
pattern = Nmer2D(n=6, d=0.2)
println("  Created pattern: ", typeof(pattern))

# Test 4: Create a molecule
println("Test 4: Testing molecule creation...")
molecule = GenericFluor(γ=10000.0, q=[0 10; 1e-2 0])
println("  Created molecule: ", typeof(molecule))

# Test 5: Run a small simulation
println("Test 5: Running a small simulation...")
try
    smld_true, smld_model, smld_noisy = simulate(
        params,
        pattern=pattern,
        molecule=molecule,
        camera=camera,
        ndatasets=1,
        nframes=10  # Small number for quick test
    )
    println("  Simulation successful!")
    println("  - True emitters: ", length(smld_true.emitters))
    println("  - Model emitters: ", length(smld_model.emitters))
    println("  - Noisy emitters: ", length(smld_noisy.emitters))
catch e
    println("  Simulation failed with error: ", e)
end

# Test 6: Test PSF integration
println("Test 6: Testing PSF integration with MicroscopePSFs...")
try
    # Create a simple emitter
    emitter = Emitter2D(1.0, 2.0, 1000.0, 1, 1)
    
    # Create a simple PSF implementation that just returns a Gaussian function
    # We'll define a minimal AbstractPSF implementation here
    struct SimplePSF <: MicroscopePSFs.AbstractPSF
        σ::Float64
    end
    
    # Define evaluate method for our simple PSF
    function MicroscopePSFs.evaluate(psf::SimplePSF, x, y, z=0.0)
        σ2 = psf.σ^2
        return exp(-(x^2 + y^2) / (2 * σ2)) / (2π * σ2)
    end
    
    # Create a PSF
    psf = SimplePSF(0.15)  # 150nm PSF width
    
    # Generate a small image
    img = gen_image(SMLD([emitter], 1, 1, 1, camera), psf, 1)
    
    println("  Successfully created image with dimensions: ", size(img))
catch e
    println("  PSF integration failed with error: ", e)
end

println("All tests completed!")
