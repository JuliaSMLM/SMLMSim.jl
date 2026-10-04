# Stepper throughput on the MicroscopeAdapt configuration (design section 7a): one `step!` plus one
# `scmos_noise!`, one thread, 256 x 256 px, about 250 in-focus emitters with blinking, bleaching, births and
# dimers, a focused-spot excitation with a switching spot, n_sub = 8, and the Cell9-like background (structured
# pattern, about 20 OOF blobs at sigma 5 px, about 30 haze blobs at sigma 2.5 px). Section 7a's render budget has
# every in-focus emitter and blob rendering in every sub-step, so the in-focus and blob populations are on 98% of
# the time (q = [-1 1; 50 -50]: still blinking, duty 50/51) and the spot switches on and off inside the measured
# frames. Prints the median time per stage and the allocations. The acceptance is a median of at most 5 ms for
# step! plus the noise call, measured on kitt.
#
# The stage costs are marginal: each is the median of the full world minus the median of the same world with that
# part left out (no dimers, no background, no blob and haze populations, UniformExcitation), not a timer inside
# `step!`. The emitting-row counts are the mean per frame of :mol and :oof rows with photons > 0.
#
# Usage: julia --project=dev dev/benchmark_stepper.jl
#        (test/test_closedloop.jl includes this file: closedloop/zero_alloc_full runs bench_setup, and with
#        SMLMSIM_BENCH=1 closedloop/benchmark asserts the bound and the workload)

using SMLMSim, MicroscopePSFs, Random, Statistics

const BENCH_T = 0.01          # exposure, s
const BENCH_PX = 0.1          # pixel size, um
const BENCH_N = 256

function bench_setup(seed::Integer = 1; dimers::Bool = true, background::Bool = true, oof::Bool = true)
    camera = SCMOSCamera(BENCH_N, BENCH_N, BENCH_PX, 1.6; offset = 100.0, gain = 2.0, qe = 1.0)
    L = BENCH_N * BENCH_PX
    q = [-1.0 1.0; 50.0 -50.0]
    mol = Population(name = :mol, density = 250 / L^2, lifetime = 3.0,
        mobility = [(0.67, 0.37), (0.33, 0.075)], fluor = GenericFluor(; γ = 9500.0, q = q), budget = 2e5,
        psf = GaussianPSF(0.0936))
    blobs = Population(name = :oof, layer = :oof, density = 20 / L^2, lifetime = 0.1,
        mobility = [(1.0, 0.25)], brightness_sigma = 0.3, fluor = GenericFluor(; γ = 5000.0, q = q),
        budget = 2e4, z = (0.5, 1.0), psf = GaussianPSF(5 * BENCH_PX), binds = false)
    haze = Population(name = :haze, layer = :oof, density = 30 / L^2, lifetime = 0.002,
        fluor = GenericFluor(; γ = 25000.0, q = zeros(1, 1)), z = (0.5, 1.0),
        psf = GaussianPSF(2.5 * BENCH_PX), binds = false)
    # every world keeps the full world's box, so a reduced world holds the same signal molecules
    margin = SMLMSim.Stepper.default_margin(Population[mol, blobs, haze])
    pops = oof ? Population[mol, blobs, haze] : Population[mol]
    bg = background ? BackgroundModel(level = 87.7, jitter = 0.025, contrast = 0.1, feature_size = 0.8,
                                      correlation_time = 0.3, illumination_width = 30.0) : nothing
    dk = dimers ? DimerKinetics(k_on = 50.0, r_react = 0.05, k_off = 0.5, D_rot = 1.0, d_dimer = 0.02) : nothing
    world = SimWorld(Random.Xoshiro(seed), camera, pops; n_sub = 8, background = bg, dimers = dk,
                     merge_radius = 0.25, margin)
    excitation = SpotExcitation(; base = 1.0, spots = [
        Spot(; x = 0.3L, y = 0.3L, σ = 0.0934, gain = 5.0, z_R = 0.23),
        Spot(; x = 0.7L, y = 0.5L, σ = 0.0934, gain = 5.0, z_R = 0.23, t_on = 2.013, t_off = 2.467)])
    return world, camera, excitation
end

# the :mol and :oof (blob) rows of the last frame that emitted photons
function bench_emitting(world)
    nmol = 0
    noof = 0
    for r in SMLMSim.frame_truth(world)
        r.photons > 0 || continue
        name = world.pops[r.pop].p.name
        name === :mol && (nmol += 1)
        name === :oof && (noof += 1)
    end
    return nmol, noof
end

# Median seconds of `step!` per frame (and of `scmos_noise!` when `noise`) over `nframes` frames after `nwarm`
# warm-up frames, for the world `bench_setup(seed; kw...)`; the warm-up of 200 frames puts the measured frames at
# 2.01 s on, so they include both spot switches (2.013 s and 2.467 s)
function bench_time(seed::Integer, exc_kind::Symbol; nwarm::Int, nframes::Int, kw...)
    world, camera, excitation = bench_setup(seed; kw...)
    exc = exc_kind === :spot ? excitation : SMLMSim.UniformExcitation()
    noise_rng = Random.Xoshiro(seed + 1)
    dst = zeros(BENCH_N, BENCH_N)
    t_step = Float64[]; t_noise = Float64[]
    for k in 1:nwarm
        SMLMSim.step!(world, (k - 1) * BENCH_T, k * BENCH_T, exc)
        scmos_noise!(noise_rng, copyto!(dst, world.expected), camera)
    end
    for k in nwarm+1:nwarm+nframes
        t0 = time_ns(); SMLMSim.step!(world, (k - 1) * BENCH_T, k * BENCH_T, exc); t1 = time_ns()
        push!(t_step, (t1 - t0) * 1e-9)
        t0 = time_ns(); scmos_noise!(noise_rng, copyto!(dst, world.expected), camera); t1 = time_ns()
        push!(t_noise, (t1 - t0) * 1e-9)
    end
    return t_step, t_noise
end

# Medians per stage, the allocations over 100 frames (each one `step!`, a `frame_truth` sweep and `scmos_noise!`)
# and the mean emitting rows per frame, for the full world; the same timing loop for the reduced worlds
function bench_run(; nwarm::Int = 200, nframes::Int = 300, seed::Integer = 1)
    world, camera, excitation = bench_setup(seed)
    noise_rng = Random.Xoshiro(seed + 1)
    dst = zeros(BENCH_N, BENCH_N)
    t_step = Float64[]; t_noise = Float64[]; t_uni = Float64[]
    a_step = Int[]; a_sweep = Int[]; a_noise = Int[]
    nmol = 0; noof = 0
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
    bench_emitting(world)
    for k in nwarm+2nframes+1:nwarm+2nframes+100
        push!(a_step, @allocated frame(k, excitation))
        push!(a_sweep, @allocated bench_emitting(world))
        push!(a_noise, @allocated scmos_noise!(noise_rng, copyto!(dst, world.expected), camera))
        m, o = bench_emitting(world)
        nmol += m; noof += o
    end
    reduced(; kw...) = median(bench_time(seed, :spot; nwarm, nframes, kw...)[1])
    both = t_step .+ t_noise
    return (; step = median(t_step), noise = median(t_noise), both = median(both), uniform_step = median(t_uni),
            no_dimers = reduced(; dimers = false), no_background = reduced(; background = false),
            no_oof = reduced(; oof = false),
            alloc_step = maximum(a_step), alloc_sweep = maximum(a_sweep), alloc_noise = maximum(a_noise),
            emitting_mol = nmol / 100, emitting_oof = noof / 100)
end

function bench_report(r)
    ms(x) = string(round(1e3 * x; digits = 3), " ms")
    println("stage                                   median")
    println("step! (SpotExcitation, full world)      ", ms(r.step))
    println("step! (UniformExcitation)               ", ms(r.uniform_step))
    println("step! (no dimers)                       ", ms(r.no_dimers))
    println("step! (no background)                   ", ms(r.no_background))
    println("step! (no blob and haze populations)    ", ms(r.no_oof))
    println("scmos_noise!                            ", ms(r.noise))
    println("step! + scmos_noise! (median of sum)    ", ms(r.both), "   (budget 5 ms)")
    println("marginal cost of (full minus reduced)")
    println("  dimers                                ", ms(r.step - r.no_dimers))
    println("  background                            ", ms(r.step - r.no_background))
    println("  blobs and haze                        ", ms(r.step - r.no_oof))
    println("  spot excitation                       ", ms(r.step - r.uniform_step))
    println("emitting rows per frame (photons > 0)   :mol ", round(r.emitting_mol; digits = 1),
            ", :oof blobs ", round(r.emitting_oof; digits = 1))
    println("max bytes allocated per call            step! ", r.alloc_step, ", frame_truth sweep ", r.alloc_sweep,
            ", noise ", r.alloc_noise)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("threads: ", Threads.nthreads(), "   ", Sys.cpu_info()[1].model)
    bench_report(bench_run())
end
