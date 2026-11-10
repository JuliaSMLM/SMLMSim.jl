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

# Golden values, identical on Julia 1.10.11, 1.11.9, 1.12.7 and 1.13.0 with Distributions
# 0.25.134, SMLMData 0.7.1 and MicroscopePSFs 0.5.7. BR_GOLDEN is main's output and random stream
# before the keyword existed, recorded at b016876 (the base of this branch). BR_GOLDEN_UDRAW is
# the output and stream with a factor drawn, recorded at 3a816c6 (the first commit of the
# keyword). To re-record, only when a toolchain change moves the random numbers: in a temporary
# environment that develops SMLMSim at that commit (a `git archive` of it) and adds Distributions
# and Random, evaluate this file's definitions above BR_GOLDEN, call each case as P1 or P1b does,
# and print br_golden_fingerprint(em, next) of its (em, md, next), and br_golden_factors(md) for
# P1b.
const BR_GOLDEN_FLUOR = GenericFluor(photons=1e4, k_off=50.0, k_on=5.0)

# A distribution that draws exactly one rand(rng), so a golden with it depends only on Julia's
# random stream and not on a Distributions sampler
struct UDraw <: ContinuousUnivariateDistribution end
Base.minimum(::UDraw) = 0.5
Base.maximum(::UDraw) = 1.5
Random.rand(rng::AbstractRNG, ::UDraw) = 0.5 + rand(rng)

function br_golden_simulate(; seed=20261108, kw...)
    Random.seed!(seed)
    params = StaticSMLMConfig(density=2.0, nframes=40, ndatasets=3, minphotons=50)
    noisy, info = simulate(params; pattern=Nmer2D(n=6, d=0.2), camera=IdealCamera(16, 16, 0.1),
                           molecule=BR_GOLDEN_FLUOR, kw...)
    return noisy.emitters, info.smld_model.metadata, rand()
end

function br_golden_kinetic(; kw...)
    emitters = [Emitter2DFit{Float64}(0.1 * k, 0.2 * k, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0;
                                      σ_xy=0.0, frame=1, dataset=1,
                                      track_id=isodd(k) ? k : 0, id=k) for k in 1:8]
    smld = BasicSMLD(emitters, IdealCamera(32, 32, 0.1), 1, 1, Dict{String,Any}())
    Random.seed!(20261109)
    out = kinetic_model(smld, GenericFluor(1e4, [-30.0 30.0; 5.0 -5.0]), 40, 20.0;
                        ndatasets=2, minphotons=0.0, kw...)
    return out.emitters, out.metadata, rand()
end

# The record count, every record's photons and x summed in order (foldl: no reassociation),
# the first, middle and last records, and the next rand() after the call.
function br_golden_fingerprint(em, next)
    n = length(em)
    rec(e) = (x=e.x, y=e.y, photons=e.photons, frame=e.frame, dataset=e.dataset, track_id=e.track_id)
    return (n=n, sum_photons=foldl(+, (e.photons for e in em)), sum_x=foldl(+, (e.x for e in em)),
            first=rec(em[1]), middle=rec(em[cld(n, 2)]), last=rec(em[n]), next_rand=next)
end

# The drawn factors in (dataset, track_id) order: their count, their sum and the first three
function br_golden_factors(md)
    f = md["brightness_factors"]
    ks = sort!(collect(keys(f)))
    return (n=length(ks), sum=foldl(+, (f[k] for k in ks)), first=[k => f[k] for k in ks[1:3]])
end

const BR_GOLDEN = Dict(
    "simulate" => (n = 546, sum_photons = 75244.28196763992, sum_x = 310.151949079229,
        first = (x = 0.1497072990457574, y = 1.2518301625993002, photons = 200.00000000000017, frame = 14, dataset = 1, track_id = 1),
        middle = (x = 1.278725603652781, y = 1.3590819118278106, photons = 89.54875394927898, frame = 34, dataset = 2, track_id = 17),
        last = (x = 0.074529134241171, y = 0.2073023930191488, photons = 200.00000000000017, frame = 40, dataset = 3, track_id = 38),
        next_rand = 0.3387342381109394),
    "kinetic_model" => (n = 182, sum_photons = 38070.39582538464, sum_x = 81.9999999999999,
        first = (x = 0.1, y = 0.2, photons = 49.97272395235577, frame = 3, dataset = 1, track_id = 1),
        middle = (x = 0.8, y = 1.6, photons = 500.00000000000045, frame = 31, dataset = 1, track_id = 11),
        last = (x = 0.8, y = 1.6, photons = 110.4056530080455, frame = 39, dataset = 2, track_id = 11),
        next_rand = 0.9483899313031025),
)

# With a factor drawn: simulate's seed is 20261110, as the default seed's noise draws differ on
# Julia 1.10 at a tail value
const BR_UDRAW_SEED = Dict("simulate" => (seed=20261110,), "kinetic_model" => NamedTuple())
const BR_GOLDEN_UDRAW = Dict(
    "simulate" => (
        fingerprint = (n = 708, sum_photons = 100120.29923255455, sum_x = 513.2092957062177,
            first = (x = 1.3040029682264414, y = 0.9778276114159444, photons = 176.85056688188715, frame = 11, dataset = 1, track_id = 1),
            middle = (x = 0.412539702938199, y = 0.5230549062482032, photons = 198.80528802814212, frame = 1, dataset = 2, track_id = 20),
            last = (x = 0.14485560368013398, y = 1.0135446628943747, photons = 176.38217937275098, frame = 25, dataset = 3, track_id = 48),
            next_rand = 0.15086413498368334),
        factors = (n = 144, sum = 139.00315949868872,
            first = [(1, 1) => 1.4929194575475249, (1, 2) => 0.6015207615004396, (1, 3) => 1.3445548041903552])),
    "kinetic_model" => (
        fingerprint = (n = 216, sum_photons = 43895.18309388006, sum_x = 99.59999999999994,
            first = (x = 0.1, y = 0.2, photons = 115.49340469378866, frame = 24, dataset = 1, track_id = 1),
            middle = (x = 0.8, y = 1.6, photons = 113.22776038542882, frame = 17, dataset = 1, track_id = 11),
            last = (x = 0.8, y = 1.6, photons = 212.4894930518612, frame = 40, dataset = 2, track_id = 11),
            next_rand = 0.47901542052010293),
        factors = (n = 16, sum = 15.546548980094808,
            first = [(1, 1) => 0.9353995372199828, (1, 3) => 1.3914885169109128, (1, 5) => 0.9785612024224327])),
)

same_records(a, b) = length(a) == length(b) &&
    all(fn -> isequal(getfield.(a, fn), getfield.(b, fn)), fieldnames(eltype(a)))

# A distribution that reports minimum 0 but draws below it
struct NegativeDraw <: ContinuousUnivariateDistribution end
Base.minimum(::NegativeDraw) = 0.0
Base.maximum(::NegativeDraw) = Inf
Random.rand(::AbstractRNG, ::NegativeDraw) = -0.5

@testset "Molecule brightness" begin
    @testset "P1: absent keyword gives main's output and random stream ($name)" for (name, case) in
            (("simulate", br_golden_simulate), ("kinetic_model", br_golden_kinetic))
        ref = BR_GOLDEN[name]
        em0 = nothing
        for kw in (NamedTuple(), (brightness=nothing,))
            em, md, next = case(; kw...)
            @test br_golden_fingerprint(em, next) == ref
            @test !haskey(md, "brightness")
            @test !haskey(md, "brightness_factors")
            em0 = em
        end
        # rand(rng, ::Dirac) returns its value and draws no random number (Distributions
        # src/univariate/discrete/dirac.jl), so with these Float64 models a unit factor gives
        # main's records and stream.
        em, md, next = case(; brightness=Dirac(1.0))
        @test br_golden_fingerprint(em, next) == ref
        @test same_records(em, em0)
        @test haskey(md, "brightness_factors")
    end

    @testset "P1b: a drawn factor: draw order and random stream ($name)" for (name, case) in
            (("simulate", br_golden_simulate), ("kinetic_model", br_golden_kinetic))
        # One rand() per molecule and dataset, just before that molecule's initial state: a
        # factor drawn anywhere else changes the records and the stream
        em, md, next = case(; BR_UDRAW_SEED[name]..., brightness=UDraw())
        @test br_golden_fingerprint(em, next) == BR_GOLDEN_UDRAW[name].fingerprint
        @test br_golden_factors(md) == BR_GOLDEN_UDRAW[name].factors
        @test md["brightness"] == UDraw()
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

    @testset "P3d: a later run without brightness drops earlier factors, not a caller's values" begin
        bright_run = br_run(br_start(10), Dirac(2.0))
        second = kinetic_model(bright_run, BR_FLUOR, 1, BR_FRAMERATE; minphotons=0.0, state1=1)
        @test !haskey(second.metadata, "brightness")
        @test !haskey(second.metadata, "brightness_factors")
        own = br_start(10)
        own.metadata["brightness"] = "a caller's note"
        own.metadata["brightness_factors"] = [1.0, 2.0]
        third = kinetic_model(own, BR_FLUOR, 1, BR_FRAMERATE; minphotons=0.0, state1=1)
        @test third.metadata["brightness"] == "a caller's note"
        @test third.metadata["brightness_factors"] == [1.0, 2.0]
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

    @testset "P6b: every drawn factor is checked ($(repr(bad)))" for bad in (NaN, Inf)
        # Minimum 1.0, so the up-front check passes. Each of 64 molecules draws the bad value
        # with probability 1/2: a run that never draws it has probability 2^-64.
        d = DiscreteNonParametric([1.0, bad], [0.5, 0.5])
        @test minimum(d) == 1.0
        Random.seed!(12)
        err = try br_run(br_start(64), d); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("drew the factor $(repr(bad))", err.msg)
        @test occursin("DiscreteNonParametric", err.msg)
        Random.seed!(12)
        err = try
            simulate(StaticSMLMConfig(nframes=1, minphotons=0); starting_conditions=br_start(64),
                     camera=IdealCamera(32, 32, 0.1), molecule=BR_FLUOR, state1=1, brightness=d)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("drew the factor $(repr(bad))", err.msg)
        @test occursin("DiscreteNonParametric", err.msg)
    end

    @testset "P6c: a negative draw is rejected even when the minimum is not" begin
        err = try br_run(br_start(3), NegativeDraw()); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("drew the factor -0.5", err.msg)
    end

    @testset "P6d: a factor whose photons leave Float64 is rejected" begin
        # Checked on the molecule's trace: a photon count that overflows (1e305 × 1e4, also for a
        # BigFloat model; a frame time that rounds a full frame past the largest Float64; a later
        # frame). A rate that underflows is a molecule that is in effect dark, and is accepted.
        on = zeros(2, 2)  # a molecule that stays in state 1
        rejected(f, b; nframes=1, framerate=BR_FRAMERATE) = try
            kinetic_model(br_start(1), f, nframes, framerate; state1=1, minphotons=0.0, brightness=Dirac(b))
            ""
        catch e
            e isa ArgumentError ? e.msg : "not an ArgumentError: $e"
        end
        in_trace = "not finite as a Float64"
        @test occursin(in_trace, rejected(GenericFluor(1e4, on), 1e305))
        @test occursin(in_trace, rejected(GenericFluor(big"1e4", on), 1e305; framerate=1.0))
        @test occursin(in_trace, rejected(GenericFluor(1e4, on), 1.6179238213760842e304; framerate=0.9))
        @test occursin(in_trace, rejected(GenericFluor(1e4, on), 5.3930794045869466e303; nframes=4, framerate=0.3))
        # the message names the model's rate too, so an invalid model shows as one
        msg = rejected(GenericFluor(NaN, on), 1.0)
        @test occursin("drew the factor 1.0", msg) && occursin("γ = NaN", msg)
        @test rejected(GenericFluor(1e4, on), 0.0) == ""  # b = 0: a rate of 0 is allowed
        @test rejected(GenericFluor(1e-200, on), 1e-200) == ""  # the rate underflows: dark, accepted
        # the same through static simulate
        sim_rejected(b; framerate=50.0) = try
            Random.seed!(3)
            simulate(StaticSMLMConfig(nframes=1, minphotons=0, framerate=framerate); starting_conditions=br_start(1),
                     camera=IdealCamera(32, 32, 0.1), molecule=GenericFluor(1e4, on), state1=1, brightness=Dirac(b))
            ""
        catch e
            e isa ArgumentError ? e.msg : "not an ArgumentError: $e"
        end
        @test occursin(in_trace, sim_rejected(1e305))
        @test occursin(in_trace, sim_rejected(1.6179238213760842e304; framerate=0.9))
        @test length(br_run(br_start(1), Dirac(1e300)).emitters) == 1  # a large factor that fits is accepted
    end

    @testset "P6e: a large distribution is cut short in the message ($(repr(bad)))" for bad in (NaN, 1e305)
        # 5001 support points print as about 150,000 characters; NaN is caught when drawn, 1e305
        # on the trace
        d = DiscreteNonParametric([collect(1.0:5000.0); bad], [fill(0.5 / 5000, 5000); 0.5])
        Random.seed!(12)
        err = try br_run(br_start(64), d); nothing catch e e end
        @test err isa ArgumentError
        @test occursin("drew the factor $(repr(bad))", err.msg)
        @test occursin("DiscreteNonParametric", err.msg) && occursin("characters in all", err.msg)
        @test length(err.msg) < 600
    end
end
