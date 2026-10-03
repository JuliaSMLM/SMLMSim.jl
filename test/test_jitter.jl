using SMLMSim, Test, Statistics, Random, StableRNGs, MicroscopePSFs, Distributions
using SMLMSim.Stepper: _add_emitter!, _advance!

# 32x32 pixel camera, 0.1 um pixels: FOV 0.0..3.2 um
cam32() = IdealCamera(1:32, 1:32, 0.1)
one_state(γ) = GenericFluor(; γ, q=zeros(1, 1))
two_state(γ, k_off, k_on) = GenericFluor(; γ, q=[-k_off k_off; k_on -k_on])
jpop(; kw...) = Population(; density=100.0, fluor=one_state(8200.0), psf=GaussianPSF(0.1), kw...)

# Kolmogorov-Smirnov statistic of v against distribution d
function ks_stat(v, d)
    x = sort(v)
    n = length(x)
    return maximum(max(i / n - cdf(d, x[i]), cdf(d, x[i]) - (i - 1) / n) for i in 1:n)
end

# log-brightness multiplier of every present emitter, by id
lj_by_id(w) = (ps = w.pops[1]; Dict(ps.id[i] => ps.lj[i] for i in 1:ps.n))

# pooled lag-1 correlation of the multiplier over consecutive exposures of the same id
function lag1_cor(w, starts, T)
    a, b = Float64[], Float64[]
    prev = nothing
    for t in starts
        SMLMSim.step!(w, t, t + T)
        cur = lj_by_id(w)
        if prev !== nothing
            for (id, x) in cur
                haskey(prev, id) && (push!(a, prev[id]); push!(b, x))
            end
        end
        prev = cur
    end
    return cor(a, b), length(a), prev
end

function step_allocs(w, k0, k1, exc::E) where {E}
    total = 0
    for k in k0:k1
        total += @allocated SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01, exc)
    end
    return total
end

@testset "stepper/jitter_keywords" begin
    p = jpop()
    @test p.brightness_jitter == 0.0 && p.jitter_time == 0.01
    for kw in ((brightness_jitter=-0.1,), (brightness_jitter=Inf,), (jitter_time=0.0,), (jitter_time=-1.0,))
        @test_throws ArgumentError jpop(; kw...)
    end
    @test jpop(; jitter_time=Inf).jitter_time == Inf
    # brightness_jitter = 0 draws nothing: identical to a population built without the keywords
    mk(; kw...) = SimWorld(StableRNG(7), cam32(),
                           [jpop(; fluor=two_state(1e4, 20.0, 20.0), lifetime=0.1, mobility=[(1.0, 0.1)], kw...)]; n_sub=4)
    w1, w2 = mk(), mk(; brightness_jitter=0.0, jitter_time=0.5)
    n_same = 0
    for k in 1:50
        SMLMSim.step!(w1, (k - 1) * 0.01, k * 0.01)
        SMLMSim.step!(w2, (k - 1) * 0.01, k * 0.01)
        n_same += w1.expected == w2.expected
    end
    @test n_same == 50
end

@testset "stepper/jitter_emission" begin
    # the emission rate of emitter i is γ_i exp(X_i): one one-state emitter over a full sub-step
    w = SimWorld(StableRNG(3), cam32(), [jpop(; density=0.0, brightness_jitter=0.3)]; n_sub=1)
    ps = w.pops[1]
    i = _add_emitter!(w, ps, 0.0)
    ps.lj[i] = 0.25
    e, alive = _advance!(w, ps, i, 0.0, 0.01, 0.0, UniformExcitation())
    @test alive
    @test e ≈ 8200.0 * exp(0.25) * 0.01 rtol = 1e-12
end

@testset "stepper/jitter_stats" begin
    s, τ, T = 0.3, 0.01, 0.01
    # (i) stationary at birth and after 20 back-to-back exposures; (ii) lag-1 correlation exp(-T/τ)
    w = SimWorld(StableRNG(101), cam32(), [jpop(; brightness_jitter=s, jitter_time=τ)]; n_sub=1)
    x0 = collect(values(lj_by_id(w)))
    @test length(x0) > 500
    @test ks_stat(x0, Normal(0, s)) < 1.628 / sqrt(length(x0))
    r, P, last = lag1_cor(w, [(k - 1) * T for k in 1:21], T)
    @test abs(r - exp(-T / τ)) < 3 / sqrt(P)
    xl = collect(values(last))
    @test ks_stat(xl, Normal(0, s)) < 1.628 / sqrt(length(xl))
    # (iii) a 0.01 s gap before each exposure: the gap sub-steps advance X too, lag-1 exp(-2T/τ)
    w = SimWorld(StableRNG(102), cam32(), [jpop(; brightness_jitter=s, jitter_time=τ)]; n_sub=1)
    r, P, _ = lag1_cor(w, [0.02 * (k - 1) for k in 1:21], T)
    @test abs(r - exp(-2T / τ)) < 3 / sqrt(P)
    # (iv) births and departures (swap-remove carries X): lag-1 still exp(-T/τ)
    w = SimWorld(StableRNG(104), cam32(), [jpop(; brightness_jitter=s, jitter_time=τ, lifetime=0.05)]; n_sub=1)
    r, P, _ = lag1_cor(w, [(k - 1) * T for k in 1:21], T)
    @test abs(r - exp(-T / τ)) < 3 / sqrt(P)
end

@testset "stepper/jitter_zero_alloc" begin
    w = SimWorld(StableRNG(5), cam32(), [jpop(; brightness_jitter=0.3, jitter_time=0.02)]; n_sub=4)
    step_allocs(w, 1, 20, UniformExcitation())
    @test step_allocs(w, 21, 120, UniformExcitation()) == 0
    ev = EvanescentExcitation(; depth=0.1, stray=0.05)
    step_allocs(w, 121, 140, ev)
    @test step_allocs(w, 141, 240, ev) == 0
end

@testset "stepper/evanescent" begin
    e = EvanescentExcitation(; depth=0.1, stray=0.1)
    @test e(0.0, 0.0, 0.0, 0.0) == 1.0
    @test e(0.0, 0.0, -0.3, 0.0) == 1.0
    @test e(1.0, 2.0, 0.2, 5.0) ≈ 0.1 + 0.9 * exp(-2) rtol = 1e-14
    d = EvanescentExcitation()
    @test d.depth == 0.1 && d.stray == 0.0
    @test_throws ArgumentError EvanescentExcitation(; depth=0.0)
    @test_throws ArgumentError EvanescentExcitation(; depth=-1.0)
    @test_throws ArgumentError EvanescentExcitation(; stray=1.5)
    @test_throws ArgumentError EvanescentExcitation(; stray=-0.1)
    @test next_switch(e, 0.3) == Inf
    # one-state emitters at z = 0.2: the image is exp(-2) times the uniform one (the same draws)
    mk() = SimWorld(StableRNG(5), cam32(), [jpop(; density=50.0, z=(0.2, 0.2))]; n_sub=2)
    wa, wb = mk(), mk()
    SMLMSim.step!(wa, 0.0, 0.01, UniformExcitation())
    SMLMSim.step!(wb, 0.0, 0.01, EvanescentExcitation(; depth=0.1))
    @test sum(wa.expected) > 0
    @test sum(wb.expected) ≈ exp(-2) * sum(wa.expected) rtol = 1e-12
end
