using SMLMSim, Test, Statistics, Random, Distributions

# one record per frame (exposure = frame period = dt), 1000 photons per record without jitter
jcfg(; kw...) = DiffusionSMLMConfig(; density=40.0, box_size=2.0, dt=0.01, t_max=0.21,
                                    camera_framerate=100.0, camera_exposure=0.01, kw...)

function ks_stat(v, d)
    x = sort(v)
    n = length(x)
    return maximum(max(i / n - cdf(d, x[i]), cdf(d, x[i]) - (i - 1) / n) for i in 1:n)
end

logX(e) = log(e.photons / 1000.0)

@testset "diffusion/jitter_keywords" begin
    p = DiffusionSMLMConfig()
    @test p.brightness_jitter == 0.0 && p.jitter_time == 0.01
    @test p.z_range == (0.0, 0.0) && p.excitation === nothing
    for kw in ((brightness_jitter=-0.1,), (brightness_jitter=Inf,), (jitter_time=0.0,), (jitter_time=-1.0,),
               (z_range=(0.3, 0.1),), (z_range=(-Inf, Inf),), (z_range=(0.0, NaN),), (z_range=(0.0, Inf),), (excitation=EvanescentExcitation(), ndims=3))
        @test_throws ArgumentError DiffusionSMLMConfig(; kw...)
    end
    # defaults draw nothing: bit-identical to a config with jitter 0 and no excitation
    Random.seed!(11); a, _ = simulate(jcfg())
    Random.seed!(11); b, _ = simulate(jcfg(; brightness_jitter=0.0, jitter_time=0.5, z_range=(0.1, 0.4)))
    @test a.emitters == b.emitters
end

@testset "diffusion/jitter_stats" begin
    s, τ = 0.3, 0.01
    Random.seed!(21)
    smld, _ = simulate(jcfg(; brightness_jitter=s, jitter_time=τ))
    byframe = Dict{Int,Vector{Float64}}()
    X = Dict{Tuple{Int,Int},Float64}()
    for e in smld.emitters
        push!(get!(() -> Float64[], byframe, e.frame), logX(e))
        X[(e.track_id, e.frame)] = logX(e)
    end
    # stationary Normal(0, s) in the first and last frame
    for f in (1, smld.n_frames)
        v = byframe[f]
        @test length(v) > 100
        @test ks_stat(v, Normal(0, s)) < 1.628 / sqrt(length(v))
    end
    # lag-1 (one dt) correlation exp(-dt/τ), pooled over molecules
    a = Float64[]; b = Float64[]
    for ((id, f), x) in X
        haskey(X, (id, f + 1)) && (push!(a, x); push!(b, X[(id, f + 1)]))
    end
    @test abs(cor(a, b) - exp(-0.01 / τ)) < 3 / sqrt(length(a))
end

@testset "diffusion/evanescent" begin
    e = EvanescentExcitation(; depth=0.1, stray=0.1)
    @test e(0.0, 0.0, 0.0, 0.0) == 1.0
    @test e(0.0, 0.0, -0.3, 0.0) == 1.0
    @test e(1.0, 2.0, 0.2, 5.0) ≈ 0.1 + 0.9 * exp(-2) rtol = 1e-14
    d = EvanescentExcitation()
    @test d.depth == 0.1 && d.stray == 0.0
    @test_throws ArgumentError EvanescentExcitation(; depth=0.0)
    @test_throws ArgumentError EvanescentExcitation(; stray=1.5)
    @test_throws ArgumentError EvanescentExcitation(; stray=-0.1)
    # the converted values are checked: a depth that underflows to 0.0, NaN depth, NaN stray
    @test_throws ArgumentError EvanescentExcitation(; depth=big"1e-1000")
    @test_throws ArgumentError EvanescentExcitation(; depth=NaN)
    @test_throws ArgumentError EvanescentExcitation(; stray=NaN)
    # every molecule at z = 0.2: the same trajectories, every record scaled by exp(-2)
    Random.seed!(31); a, _ = simulate(jcfg())
    Random.seed!(31); b, _ = simulate(jcfg(; z_range=(0.2, 0.2), excitation=EvanescentExcitation(; depth=0.1)))
    @test [(r.x, r.y, r.track_id, r.frame) for r in a.emitters] == [(r.x, r.y, r.track_id, r.frame) for r in b.emitters]
    @test all(i -> isapprox(b.emitters[i].photons, exp(-2) * a.emitters[i].photons; rtol=1e-12), eachindex(a.emitters))
    # z drawn once per molecule in z_range: a constant factor per molecule, within the range's bounds
    Random.seed!(32)
    c, _ = simulate(jcfg(; z_range=(0.0, 0.3), excitation=EvanescentExcitation(; depth=0.1)))
    fac = Dict{Int,Vector{Float64}}()
    for r in c.emitters
        push!(get!(() -> Float64[], fac, r.track_id), r.photons / 1000.0)
    end
    @test all(v -> maximum(v) - minimum(v) <= 1e-12 * maximum(v), values(fac))
    @test all(v -> exp(-3) - 1e-12 <= v[1] <= 1 + 1e-12, values(fac))
    @test length(unique(round(v[1]; digits=9) for v in values(fac))) > 10
end

@testset "one EvanescentExcitation and one OU step for both paths" begin
    # The diffusion path and the stepper share Core's type and coefficients
    @test SMLMSim.Stepper.EvanescentExcitation === SMLMSim.Core.EvanescentExcitation
    @test EvanescentExcitation === SMLMSim.Core.EvanescentExcitation
    @test SMLMSim.Stepper._ou_coeffs === SMLMSim.Core._ou_coeffs
    e = EvanescentExcitation(; depth=0.2, stray=0.1)
    @test e(0.0, 0.0, 0.3, 0.0) === 0.1 + 0.9 * exp(-0.3 / 0.2)
    @test e(0.0, 0.0, -1.0, 0.0) === 1.0
    @test e(0, 0, 0.3f0, 0) === 0.1 + 0.9 * exp(-Float64(0.3f0) / 0.2)
    # the exact AR(1) coefficients the stepper used inline before
    a, b = SMLMSim.Core._ou_coeffs(0.3, 0.02, 0.001)
    @test a === exp(-0.001 / 0.02) && b === 0.3 * sqrt(-expm1(-2 * 0.001 / 0.02))
end

@testset "diffusion/legacy_fixture" begin
    # SMLMSim 0.7.2's own output for a seeded default run (registered 0.7.2 in a temporary environment): the
    # record count, sums of x, y, photons and timestamp, and the next rand() after simulate. Julia 1.11
    # changed randn's tail sampler, so the draws differ per version; 1.10 values recorded on 1.10.11 and 1.11+
    # values on 1.12.7.
    cfg = DiffusionSMLMConfig(density=40.0, box_size=2.0, diff_monomer=0.2, diff_dimer=0.1, diff_dimer_rot=0.5,
        k_off=0.5, r_react=0.05, d_dimer=0.03, dt=0.01, t_max=1.0, camera_framerate=100.0, camera_exposure=0.01)
    Random.seed!(20261003)
    smld, _ = simulate(cfg)
    e = smld.emitters
    got = (length(e), sum(r.x for r in e), sum(r.y for r in e), sum(r.photons for r in e), sum(r.timestamp for r in e))
    nxt = rand()
    gold, gold_next = VERSION >= v"1.11-" ?
        ((16000, 16178.292676401472, 16051.005437711476, 1.6e7, 7920.00000000007), 0.5483578200259479) :
        ((16000, 15815.511433569145, 16923.313494147566, 1.6e7, 7920.00000000007), 0.2860410262182773)
    @test got[1] == gold[1]
    @test all(i -> isapprox(got[i], gold[i]; rtol=1e-12), 2:5)
    @test nxt == gold_next
    @test !haskey(smld.metadata, "base_photons")
end

@testset "diffusion/ou_unrecorded_steps" begin
    # dt = 1 ms, 5 ms exposure, 10 ms frame period: X advances every dt, so the frame-start log-photon
    # correlation across one frame period is exp(-0.01/τ), not exp(-0.005/τ)
    s, τ = 0.3, 0.02
    cfg = DiffusionSMLMConfig(; density=40.0, box_size=2.0, dt=0.001, t_max=0.5, camera_framerate=100.0,
                              camera_exposure=0.005, brightness_jitter=s, jitter_time=τ)
    Random.seed!(41)
    smld, _ = simulate(cfg)
    first = Dict{Tuple{Int,Int},Tuple{Float64,Float64}}()
    for e in smld.emitters
        k = (e.track_id, e.frame)
        (!haskey(first, k) || e.timestamp < first[k][1]) && (first[k] = (e.timestamp, log(e.photons / 100.0)))
    end
    a = Float64[]; b = Float64[]
    for ((id, f), (_, x)) in first
        haskey(first, (id, f + 1)) && (push!(a, x); push!(b, first[(id, f + 1)][2]))
    end
    @test length(a) > 5000
    @test abs(cor(a, b) - exp(-0.01 / τ)) < 4 / sqrt(length(a))
    @test abs(cor(a, b) - exp(-0.005 / τ)) > 0.1

    # the two members of a dimer carry independent factors
    mk(id, partner) = SMLMSim.DiffusingEmitter2D{Float64}(1.0, 1.0, 1000.0, 0.0, 1, 1, id, :dimer, partner)
    pair = SMLMSim.DiffusingEmitter2D{Float64}[mk(1, 2), mk(2, 1)]
    pcfg = DiffusionSMLMConfig(; box_size=2.0, diff_monomer=0.0, diff_dimer=0.0, diff_dimer_rot=0.0, k_off=0.0,
                               dt=0.01, t_max=3.0, camera_framerate=100.0, camera_exposure=0.01,
                               brightness_jitter=0.3, jitter_time=0.01)
    Random.seed!(42)
    ps, _ = simulate(pcfg; starting_conditions=pair)
    @test all(r -> r.state == :dimer, ps.emitters)
    x1 = [log(r.photons / 1000.0) for r in ps.emitters if r.track_id == 1]
    x2 = [log(r.photons / 1000.0) for r in ps.emitters if r.track_id == 2]
    @test length(x1) == length(x2) > 200
    @test x1 != x2
    @test abs(cor(x1, x2)) < 0.25
end

@testset "diffusion/base_photons continuation" begin
    # Codex's case: γ = 1e5, one step per frame, no motion, no jitter, every molecule at z = 0.2, depth 0.1
    # (each record carries 1000·exp(-2) = 135.335)
    keep1(smld) = BasicSMLD([e for e in smld.emitters if e.track_id == 1], smld.camera, smld.n_frames,
                            smld.n_datasets, copy(smld.metadata))
    cfg = DiffusionSMLMConfig(; box_size=2.0, diff_monomer=0.0, dt=0.01, t_max=0.03, camera_framerate=100.0,
                              camera_exposure=0.01, z_range=(0.2, 0.2), excitation=EvanescentExcitation(; depth=0.1))
    Random.seed!(51)
    smld, _ = simulate(cfg; γ=1e5, override_count=2)
    @test all(e -> isapprox(e.photons, 1000 * exp(-2); rtol=1e-12), smld.emitters)
    @test smld.metadata["base_photons"] == Dict{Int,Float64}(1 => 1000.0, 2 => 1000.0)
    full, _ = simulate(cfg; starting_conditions=smld)
    sub = keep1(smld)
    cont, _ = simulate(cfg; starting_conditions=sub)
    @test !isempty(cont.emitters) && all(e -> e.track_id == 1, cont.emitters)
    @test all(e -> isapprox(e.photons, 1000 * exp(-2); rtol=1e-12), full.emitters)
    @test all(e -> isapprox(e.photons, 1000 * exp(-2); rtol=1e-12), cont.emitters)
    # extract_end_state copies it, unknown provenance drops it
    @test extract_end_state(sub).metadata["base_photons"] == smld.metadata["base_photons"]
    tied = BasicSMLD(vcat(smld.emitters, smld.emitters), smld.camera, smld.n_frames, smld.n_datasets, copy(smld.metadata))
    @test !haskey(extract_end_state(tied).metadata, "base_photons")
    # no modulation: the key is absent
    plain, _ = simulate(DiffusionSMLMConfig(; box_size=2.0, t_max=0.03, camera_framerate=100.0, camera_exposure=0.01))
    @test !haskey(plain.metadata, "base_photons")

    # Jitter: the filtered continuation resumes at the saved base, not at a modulated record
    jc = DiffusionSMLMConfig(; box_size=2.0, diff_monomer=0.0, dt=0.01, t_max=0.05, camera_framerate=100.0,
                             camera_exposure=0.01, brightness_jitter=0.5)
    Random.seed!(52)
    js, _ = simulate(jc; γ=1e5, override_count=2)
    @test any(e -> !isapprox(e.photons, 1000.0; rtol=1e-3), js.emitters)
    @test js.metadata["base_photons"] == Dict{Int,Float64}(1 => 1000.0, 2 => 1000.0)
    @test [e.photons for e in extract_end_state(keep1(js)).emitters] == [1000.0]
    @test [e.photons for e in extract_end_state(js).emitters] == [1000.0, 1000.0]
end
