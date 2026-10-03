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
               (z_range=(0.3, 0.1),), (excitation=EvanescentExcitation(), ndims=3))
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
