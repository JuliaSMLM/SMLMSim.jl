using SMLMSim, Test, Distributions, Statistics, Random

# Molecules that stay on: state1=1 and an off rate so small that none leaves state 1.
const BR_γ = 1e4
const BR_FRAMERATE = 20.0
const BR_FLUOR = GenericFluor(BR_γ, [-1e-12 1e-12; 1.0 -1.0])

function br_start(N; ndims=2)
    emitters = if ndims == 2
        [Emitter2DFit{Float64}(1.0, 2.0, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0;
                               σ_xy=0.0, frame=1, dataset=1, track_id=i, id=i) for i in 1:N]
    else
        [Emitter3DFit{Float64}(1.0, 2.0, 0.5, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
                               σ_xy=0.0, frame=1, dataset=1, track_id=i, id=i) for i in 1:N]
    end
    return BasicSMLD(emitters, IdealCamera(32, 32, 0.1), 1, 1, Dict{String,Any}())
end

br_run(smld, b; nframes=1, ndatasets=1, minphotons=0.0) =
    kinetic_model(smld, BR_FLUOR, nframes, BR_FRAMERATE; ndatasets=ndatasets,
                  minphotons=minphotons, state1=1, brightness=b)

# Exact mean and variance of max(Normal(μ, σ), a).
function censored_moments(μ, σ, a)
    α = (a - μ) / σ
    Φ, φ = cdf(Normal(), α), pdf(Normal(), α)
    m1 = a * Φ + μ * (1 - Φ) + σ * φ
    m2 = a^2 * Φ + (μ^2 + σ^2) * (1 - Φ) + σ * (a + μ) * φ
    return m1, m2 - m1^2
end

@testset "Molecule brightness" begin
    @testset "P1: absent keyword changes nothing" begin
        cam = IdealCamera(32, 32, 0.1)
        params = StaticSMLMConfig(density=10.0, nframes=20, ndatasets=2, minphotons=10)
        run(; kw...) = (Random.seed!(7);
            simulate(params; pattern=Nmer2D(n=4, d=0.2), camera=cam, kw...))
        s0, i0 = run()
        s1, i1 = run(brightness=nothing)
        @test length(s0.emitters) == length(s1.emitters) > 0
        @test all(fieldnames(Emitter2DFit{Float64})) do fn
            all(getfield.(s0.emitters, fn) .== getfield.(s1.emitters, fn))
        end
        @test sort(collect(keys(i0.smld_model.metadata))) == sort(collect(keys(i1.smld_model.metadata)))
        @test !haskey(i0.smld_model.metadata, "brightness")
        @test !haskey(i0.smld_model.metadata, "brightness_factors")

        # rand(rng, ::Dirac) returns its value and draws no random number (Distributions
        # src/univariate/discrete/dirac.jl), so a unit factor reproduces the default exactly.
        s2, i2 = run(brightness=Dirac(1.0))
        @test length(s0.emitters) == length(s2.emitters)
        @test all(fieldnames(Emitter2DFit{Float64})) do fn
            all(getfield.(s0.emitters, fn) .== getfield.(s2.emitters, fn))
        end
        @test haskey(i2.smld_model.metadata, "brightness_factors")
    end

    @testset "P2: factors follow the distribution" begin
        N = 20_000
        μ, σ, a = 1.0, 0.3, 0.7
        cases = (
            ("censored", censored(Normal(μ, σ); lower=a), censored_moments(μ, σ, a)),
            ("truncated", truncated(Normal(μ, σ); lower=a),
             (mean(truncated(Normal(μ, σ); lower=a)), var(truncated(Normal(μ, σ); lower=a)))),
        )
        for (name, dist, (m, v)) in cases
            @testset "$name" begin
                Random.seed!(11)
                out = br_run(br_start(N), dist)
                @test length(out.emitters) == N
                x = [e.photons * BR_FRAMERATE / BR_γ for e in out.emitters]
                se_m = sqrt(v / N)
                @test abs(mean(x) - m) < 4 * se_m
                m4 = mean((x .- mean(x)) .^ 4)
                se_v = sqrt((m4 - var(x)^2) / N)
                @test abs(var(x) - v) < 4 * se_v
                @test all(x .>= a)
                facs = out.metadata["brightness_factors"]
                @test facs isa Dict{Tuple{Int,Int},Float64}
                @test length(facs) == N
                @test out.metadata["brightness"] == dist
                @test all(e -> isapprox(e.photons,
                    facs[(e.dataset, e.track_id)] * BR_γ / BR_FRAMERATE; rtol=1e-12), out.emitters)
                if name == "censored"
                    p0 = cdf(Normal(μ, σ), a)
                    frac = count(==(a), values(facs)) / N
                    @test abs(frac - p0) < 4 * sqrt(p0 * (1 - p0) / N)
                else
                    @test count(==(a), values(facs)) == 0
                end
            end
        end
    end

    @testset "P3: one factor per molecule, fresh per dataset" begin
        N = 500
        Random.seed!(3)
        out = br_run(br_start(N), truncated(Normal(1, 0.3); lower=0.7); nframes=3, ndatasets=2)
        @test length(out.emitters) == 3 * N * 2
        same_within = all(Iterators.product(1:2, 1:N)) do (dd, id)
            ph = [e.photons for e in out.emitters if e.dataset == dd && e.track_id == id]
            length(ph) == 3 && all(isapprox.(ph, ph[1]; rtol=1e-12))
        end
        @test same_within
        facs = out.metadata["brightness_factors"]
        @test length(facs) == 2N
        # a continuous distribution: every molecule's two factors differ
        @test count(id -> facs[(1, id)] != facs[(2, id)], 1:N) == N
    end

    @testset "P3b: the rate matrix stays the model's" begin
        # A blinking model: with no brightness draw (Dirac), the chain's draws are the same, so a
        # factor of 2 must give the same frames with photons exactly doubled.
        fluor = GenericFluor(1e4, [-30.0 30.0; 5.0 -5.0])
        smld = br_start(200)
        run(b) = (Random.seed!(21);
            kinetic_model(smld, fluor, 50, BR_FRAMERATE; ndatasets=2, minphotons=0.0, brightness=b))
        e1, e2 = run(nothing).emitters, run(Dirac(2.0)).emitters
        @test length(e1) == length(e2) > 0
        @test [(e.dataset, e.track_id, e.frame) for e in e1] == [(e.dataset, e.track_id, e.frame) for e in e2]
        @test all(2 .* getfield.(e1, :photons) .== getfield.(e2, :photons))
        @test any(e -> e.photons < 1e4 / BR_FRAMERATE, e1)   # some partial frames: the model blinks
    end

    @testset "P3c: factors are keyed by the output track_id" begin
        ids = [0, 0, 5, 0]
        emitters = [Emitter2DFit{Float64}(1.0, 2.0, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0;
                                          σ_xy=0.0, frame=1, dataset=1, track_id=t, id=k)
                    for (k, t) in enumerate(ids)]
        smld = BasicSMLD(emitters, IdealCamera(32, 32, 0.1), 1, 1, Dict{String,Any}())
        Random.seed!(4)
        out = br_run(smld, truncated(Normal(1, 0.3); lower=0.7); ndatasets=2)
        facs = out.metadata["brightness_factors"]
        @test sort(collect(keys(facs))) == [(d, t) for d in 1:2 for t in 5:8]
        @test all(e -> isapprox(e.photons,
            facs[(e.dataset, e.track_id)] * BR_γ / BR_FRAMERATE; rtol=1e-12), out.emitters)
    end

    @testset "P3d: a later run without brightness drops earlier factors" begin
        bright_run = br_run(br_start(10), Dirac(2.0))
        second = kinetic_model(bright_run, BR_FLUOR, 1, BR_FRAMERATE; minphotons=0.0, state1=1)
        @test !haskey(second.metadata, "brightness")
        @test !haskey(second.metadata, "brightness_factors")
    end

    @testset "P4: 2D and 3D, starting conditions and pattern path" begin
        cam = IdealCamera(32, 32, 0.1)
        b = truncated(Normal(1, 0.3); lower=0.7)
        params = StaticSMLMConfig(density=10.0, nframes=20, ndatasets=2, minphotons=10)
        for (pat, start, T) in (
            (Nmer2D(n=4, d=0.2), br_start(30), Emitter2DFit),
            (Nmer3D(n=4, d=0.2), br_start(30; ndims=3), Emitter3DFit),
        )
            Random.seed!(5)
            _, info = simulate(StaticSMLMConfig(density=10.0, nframes=20, ndatasets=2, minphotons=10,
                                                ndims=T === Emitter2DFit ? 2 : 3);
                               pattern=pat, camera=cam, brightness=b)
            em = info.smld_model.emitters
            @test !isempty(em)
            @test eltype(em) <: T
            facs = info.smld_model.metadata["brightness_factors"]
            @test all(e -> haskey(facs, (e.dataset, e.track_id)), em)
            @test info.smld_model.metadata["brightness"] == b

            _, info = simulate(params; starting_conditions=start, camera=cam,
                               molecule=BR_FLUOR, state1=1, brightness=b)
            em = info.smld_model.emitters
            @test !isempty(em)
            @test eltype(em) <: T
            facs = info.smld_model.metadata["brightness_factors"]
            @test all(e -> haskey(facs, (e.dataset, e.track_id)), em)
        end
    end

    @testset "P5: discrete family; zero factor emits nothing" begin
        N = 2000
        Random.seed!(9)
        out = br_run(br_start(N), Poisson(2))   # minphotons = 0: the rule is photons > 0
        facs = out.metadata["brightness_factors"]
        @test facs isa Dict{Tuple{Int,Int},Float64}
        @test all(v -> v == round(v) && v >= 0, values(facs))
        zeros_ids = [k[2] for (k, v) in facs if v == 0]
        @test !isempty(zeros_ids)
        @test length(out.emitters) == N - length(zeros_ids)
        @test isdisjoint(Set(e.track_id for e in out.emitters), Set(zeros_ids))
    end

    @testset "P6: argument checks" begin
        smld = br_start(5)
        @test_throws ArgumentError br_run(smld, Normal(1, 0.3))
        err = try br_run(smld, Normal(1, 0.3)) catch e e end
        @test occursin("truncated", err.msg) && occursin("censored", err.msg)
        struct OtherFluor <: Molecule end
        @test_throws ArgumentError kinetic_model(smld, OtherFluor(), 1, 10.0; brightness=Dirac(1.0))
        # the checks come before any random draw
        Random.seed!(1); r1 = rand()
        Random.seed!(1); try br_run(smld, Normal(1, 0.3)) catch end
        @test rand() == r1
        # in simulate too, before the pattern and labeling draws
        params = StaticSMLMConfig(density=10.0, nframes=5)
        Random.seed!(1); try simulate(params; camera=IdealCamera(16, 16, 0.1), brightness=Normal(1, 0.3)) catch end
        @test rand() == r1
        @test_throws ArgumentError simulate(params; camera=IdealCamera(16, 16, 0.1), brightness=Normal(1, 0.3))
        # a NaN factor is rejected
        @test_throws ArgumentError br_run(smld, Dirac(NaN))
    end
end
