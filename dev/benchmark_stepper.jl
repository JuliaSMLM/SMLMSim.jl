# Stepper throughput on the MicroscopeAdapt configuration (design section 7a): one `step!` plus one
# `scmos_noise!`, one thread, 256 x 256 px, about 250 in-focus emitters with blinking, bleaching, dimers and
# a focused-spot excitation, n_sub = 8, and the Cell9-like background (structured pattern, about 20 OOF blobs
# at sigma 5 px, about 30 haze blobs at sigma 2.5 px). Prints the median time per stage and the allocations.
# The acceptance is a median of at most 5 ms for step! plus the noise call, measured on kitt.
#
# Usage: julia --project=dev dev/benchmark_stepper.jl
#        (test/test_closedloop.jl includes this file and asserts the bound when SMLMSIM_BENCH=1)

using SMLMSim, MicroscopePSFs, Random, Statistics

const BENCH_T = 0.01          # exposure, s
const BENCH_PX = 0.1          # pixel size, um
const BENCH_N = 256

function bench_setup(seed::Integer = 1)
    camera = SCMOSCamera(BENCH_N, BENCH_N, BENCH_PX, 1.6; offset = 100.0, gain = 2.0, qe = 1.0)
    L = BENCH_N * BENCH_PX
    q = [-50.0 50.0; 1.0 -1.0]
    mol = Population(name = :mol, density = 250 / L^2, lifetime = 3.0,
        mobility = [(0.67, 0.37), (0.33, 0.075)], fluor = GenericFluor(; γ = 9500.0, q = q), budget = 5e4,
        psf = GaussianPSF(0.0936))
    blobs = Population(name = :oof, layer = :oof, density = 20 / L^2, lifetime = 0.1,
        mobility = [(1.0, 0.25)], brightness_sigma = 0.3, fluor = GenericFluor(; γ = 5000.0, q = q),
        budget = 2e4, z = (0.5, 1.0), psf = GaussianPSF(5 * BENCH_PX), binds = false)
    haze = Population(name = :haze, layer = :oof, density = 30 / L^2, lifetime = 0.002,
        fluor = GenericFluor(; γ = 25000.0, q = zeros(1, 1)), z = (0.5, 1.0),
        psf = GaussianPSF(2.5 * BENCH_PX), binds = false)
    bg = BackgroundModel(level = 87.7, jitter = 0.025, contrast = 0.1, feature_size = 0.8,
                         correlation_time = 0.3, illumination_width = 30.0)
    dimers = DimerKinetics(k_on = 50.0, r_react = 0.05, k_off = 0.5, D_rot = 1.0, d_dimer = 0.02)
    world = SimWorld(Random.Xoshiro(seed), camera, [mol, blobs, haze]; n_sub = 8, background = bg,
                     dimers = dimers, merge_radius = 0.25)
    excitation = SpotExcitation(; base = 1.0, spots = [
        Spot(; x = 0.3L, y = 0.3L, σ = 0.0934, gain = 5.0, z_R = 0.23),
        Spot(; x = 0.7L, y = 0.5L, σ = 0.0934, gain = 5.0, z_R = 0.23, t_on = 1.0, t_off = 2.0)])
    return world, camera, excitation
end

# Median seconds per stage over `nframes` frames after `nwarm` warm-up frames, and allocated bytes per frame
function bench_run(; nwarm::Int = 300, nframes::Int = 300, seed::Integer = 1)
    world, camera, excitation = bench_setup(seed)
    noise_rng = Random.Xoshiro(seed + 1)
    dst = zeros(BENCH_N, BENCH_N)
    t_step = Float64[]; t_noise = Float64[]; t_uni = Float64[]
    a_step = Int[]; a_noise = Int[]
    frame(k, exc) = SMLMSim.step!(world, (k - 1) * BENCH_T, k * BENCH_T, exc)
    for k in 1:nwarm
        frame(k, excitation)
        scmos_noise!(noise_rng, copyto!(dst, world.expected), camera)
    end
    for k in nwarm+1:nwarm+nframes
        t0 = time_ns(); frame(k, excitation); t1 = time_ns()
        push!(t_step, (t1 - t0) * 1e-9)
        t0 = time_ns(); scmos_noise!(noise_rng, copyto!(dst, world.expected), camera); t1 = time_ns()
        push!(t_noise, (t1 - t0) * 1e-9)
    end
    for k in nwarm+nframes+1:nwarm+2nframes
        t0 = time_ns(); frame(k, SMLMSim.UniformExcitation()); t1 = time_ns()
        push!(t_uni, (t1 - t0) * 1e-9)
    end
    for k in nwarm+2nframes+1:nwarm+2nframes+20
        push!(a_step, @allocated frame(k, excitation))
        push!(a_noise, @allocated scmos_noise!(noise_rng, copyto!(dst, world.expected), camera))
    end
    both = t_step .+ t_noise
    return (; step = median(t_step), noise = median(t_noise), both = median(both), uniform_step = median(t_uni),
            alloc_step = maximum(a_step), alloc_noise = maximum(a_noise),
            n_emitters = sum(r -> world.pops[r.pop].p.layer === :signal, SMLMSim.frame_truth(world)))
end

function bench_report(r)
    ms(x) = string(round(1e3 * x; digits = 3), " ms")
    println("stage                                   median")
    println("step! (SpotExcitation)                  ", ms(r.step))
    println("step! (UniformExcitation)               ", ms(r.uniform_step))
    println("scmos_noise!                            ", ms(r.noise))
    println("step! + scmos_noise! (median of sum)    ", ms(r.both), "   (budget 5 ms)")
    println("max bytes allocated per call            step! ", r.alloc_step, ", noise ", r.alloc_noise)
    println("in-focus emitters in the last frame     ", r.n_emitters)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("threads: ", Threads.nthreads(), "   ", Sys.cpu_info()[1].model)
    bench_report(bench_run())
end
