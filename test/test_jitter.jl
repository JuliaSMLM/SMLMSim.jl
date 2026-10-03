using SMLMSim, Test, Statistics, Random, StableRNGs, MicroscopePSFs, Distributions
using SMLMSim.Stepper: _add_emitter!, _advance!
using SMLMSim: frame_truth

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

# X = log(photons/(γ T)) of every full-exposure row of the last step (NaN t_depart, and NaN t_birth or a birth at
# t_a, as the initial emitters of exposure 1 have), by id: for a one-state emitter at n_sub = 1 the photons are
# γ exp(X) T
function x_by_id(w, T; γ=8200.0)
    d = Dict{Int,Float64}()
    for r in frame_truth(w)
        (isnan(r.t_birth) || r.t_birth == w.t_a) && isnan(r.t_depart) && (d[r.id] = log(r.photons / (γ * T)))
    end
    return d
end

# X of the full-exposure rows of back-to-back or gapped exposures starting at `starts`, one Dict per exposure
function x_frames(w, starts, T)
    out = Dict{Int,Float64}[]
    for t in starts
        SMLMSim.step!(w, t, t + T)
        push!(out, x_by_id(w, T))
    end
    return out
end

# pooled lag-1 correlation of X over consecutive exposures of the same id, and the pair count
function lag1_cor(frames)
    a, b = Float64[], Float64[]
    for k in 2:length(frames), (id, x) in frames[k]
        haskey(frames[k-1], id) && (push!(a, frames[k-1][id]); push!(b, x))
    end
    return cor(a, b), length(a)
end

function count_allocs(w, exc::E, nrng, dst, cam, k) where {E}
    return @allocated begin
        SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01, exc)
        SMLMSim.scmos_noise!(nrng, copyto!(dst, w.expected), cam)
    end
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
    # X is read from the photons of full-exposure rows
    s, τ, T = 0.3, 0.01, 0.01
    # (i) stationary in exposure 1 and in exposure 21; (ii) lag-1 correlation exp(-T/τ)
    w = SimWorld(StableRNG(101), cam32(), [jpop(; brightness_jitter=s, jitter_time=τ)]; n_sub=1)
    fr = x_frames(w, [(k - 1) * T for k in 1:21], T)
    for k in (1, 21)
        xs = collect(values(fr[k]))
        @test length(xs) > 500
        @test ks_stat(xs, Normal(0, s)) < 1.628 / sqrt(length(xs))
    end
    r, P = lag1_cor(fr)
    @test abs(r - exp(-T / τ)) < 3 / sqrt(P)
    # (iii) a 0.01 s gap between exposures: the gap sub-steps advance X too, lag-1 exp(-2T/τ)
    w = SimWorld(StableRNG(102), cam32(), [jpop(; brightness_jitter=s, jitter_time=τ)]; n_sub=1)
    r, P = lag1_cor(x_frames(w, [0.02 * (k - 1) for k in 1:21], T))
    @test abs(r - exp(-2T / τ)) < 3 / sqrt(P)
    # (iv) exposure averaging at n_sub = 10: the sd of log photons is s sqrt(g) (Population docstring)
    n, ρ = 10, exp(-0.1)
    g = (n + 2 * sum((n - k) * ρ^k for k in 1:n-1)) / n^2
    @test g ≈ 0.7384782790044377 rtol = 1e-14
    w = SimWorld(StableRNG(103), cam32(), [jpop(; brightness_jitter=0.1, jitter_time=τ)]; n_sub=10)
    fr = x_frames(w, [(k - 1) * T for k in 1:21], T)
    Y = reduce(vcat, [collect(values(f)) for f in fr])
    @test abs(std(Y) / (0.1 * sqrt(g)) - 1) < 0.03
    # (v) births and departures (swap-remove carries X): lag-1 over full-exposure rows still exp(-T/τ)
    w = SimWorld(StableRNG(104), cam32(), [jpop(; brightness_jitter=s, jitter_time=τ, lifetime=0.05)]; n_sub=1)
    r, P = lag1_cor(x_frames(w, [(k - 1) * T for k in 1:21], T))
    @test abs(r - exp(-T / τ)) < 3 / sqrt(P)
end

@testset "stepper/jitter_page" begin
    # calibration to the measured within-track sd of log photons per 10 ms frame: 0.37 (Cell9) and 0.48 (Cell1),
    # baseline 0.24 without jitter (~/julia_shared_dev/LidkeLab/PPIDetect/dev/output/t15/sim_vs_real.md).
    # At n_sub = 1 and jitter_time = 0.002 the exposures are nearly independent, so the per-frame sd is brightness_jitter.
    T = 0.01
    for (target, seed) in ((0.37, 105), (0.48, 106))
        s = sqrt(target^2 - 0.24^2)
        w = SimWorld(StableRNG(seed), cam32(), [jpop(; brightness_jitter=s, jitter_time=0.002)]; n_sub=1)
        fr = x_frames(w, [(k - 1) * T for k in 1:20], T)
        J = median([std([f[id] for f in fr]) for id in keys(fr[1]) if all(f -> haskey(f, id), fr)])
        @test abs(sqrt(J^2 + 0.24^2) - target) < 0.02
    end
end

@testset "stepper/jitter_newborn" begin
    # a newborn is stationary at birth and its multiplier advances from its birth time (exact OU)
    s, τ, T = 0.3, 0.01, 0.01
    pop = Population(density=0.0, lifetime=1.0, birth_rate=2000.0, fluor=one_state(8200.0), psf=GaussianPSF(0.1),
                     brightness_jitter=s, jitter_time=τ)
    w = SimWorld(StableRNG(201), cam32(), [pop]; n_sub=1)
    Xb = Float64[]
    Z = Float64[]
    born = Dict{Int,Tuple{Float64,Float64}}()        # id => (X_b, a) of the previous exposure's births
    for k in 1:21
        t_a = (k - 1) * T
        SMLMSim.step!(w, t_a, t_a + T)
        nextborn = Dict{Int,Tuple{Float64,Float64}}()
        for r in frame_truth(w)
            if isnan(r.t_birth) && isnan(r.t_depart) && haskey(born, r.id)
                Xb0, a = born[r.id]
                push!(Z, (log(r.photons / (8200.0 * T)) - a * Xb0) / (s * sqrt(1 - a^2)))
            elseif !isnan(r.t_birth) && isnan(r.t_depart) && t_a + T - r.t_birth > 1e-6
                d = t_a + T - r.t_birth
                x = log(r.photons / (8200.0 * d))
                push!(Xb, x)
                nextborn[r.id] = (x, exp(-d / τ))
            end
        end
        born = nextborn
    end
    M = length(Z)
    @test M > 3000
    @test abs(mean(Z)) < 4 / sqrt(M)
    @test abs(var(Z) - 1) < 0.1
    @test ks_stat(Xb, Normal(0, s)) < 1.628 / sqrt(length(Xb))
end

@testset "stepper/jitter_draw_order" begin
    # the jitter pass draws after every other draw of the sub-step (motion of an :oof population here)
    p1 = Population(name=:sig, density=0.0, fluor=one_state(1000.0), psf=GaussianPSF(0.1), brightness_jitter=0.3)
    p2 = Population(name=:oof, layer=:oof, density=0.0, mobility=[(1.0, 0.1)], fluor=one_state(1000.0),
                    psf=GaussianPSF(0.39), z=(0.6, 0.6))
    w = SimWorld(StableRNG(7), cam32(), [p1, p2]; n_sub=1)
    i1 = _add_emitter!(w, w.pops[1], 0.0)
    i2 = _add_emitter!(w, w.pops[2], 0.0)
    w.pops[1].x[i1], w.pops[1].y[i1], w.pops[1].lj[i1] = 1.0, 1.0, 0.1
    w.pops[2].x[i2], w.pops[2].y[i2] = 1.6, 2.4
    r = copy(w.rng)
    ξx = randn(r); ξy = randn(r); ξj = randn(r)
    SMLMSim.step!(w, 0.0, 0.01)
    sd = sqrt(2 * 0.1 * 0.01)
    a, b = SMLMSim.Core._ou_coeffs(0.3, 0.01, 0.01)
    @test w.pops[2].x[i2] ≈ 1.6 + sd * ξx rtol = 1e-14
    @test w.pops[2].y[i2] ≈ 2.4 + sd * ξy rtol = 1e-14
    @test w.pops[1].lj[i1] ≈ a * 0.1 + b * ξj rtol = 1e-14
end

@testset "stepper/jitter_replay" begin
    mk() = SimWorld(StableRNG(7), cam32(),
                    [jpop(; fluor=two_state(1e4, 20.0, 20.0), lifetime=0.1, budget=5000.0, mobility=[(1.0, 0.1)],
                          brightness_jitter=0.3)]; n_sub=4)
    a, b = mk(), mk()
    ok = true
    for k in 1:50
        ea = SMLMSim.step!(a, (k - 1) * 0.01, k * 0.01)
        eb = SMLMSim.step!(b, (k - 1) * 0.01, k * 0.01)
        ok &= ea == eb && collect(frame_truth(a)) == collect(frame_truth(b))
    end
    @test ok
end

@testset "stepper/jitter_growth" begin
    # the capacity starts at 64 and grows: lj grows with the other vectors and the correlation holds
    w = SimWorld(StableRNG(110), cam32(), [jpop(; density=0.0, lifetime=0.05, birth_rate=2000.0, brightness_jitter=0.3)]; n_sub=1)
    ps = w.pops[1]
    T = 0.01
    fr = x_frames(w, [(k - 1) * T for k in 1:21], T)
    @test length(ps.lj) == length(ps.x)
    @test length(ps.x) > 64
    r, P = lag1_cor(fr[10:21])
    @test abs(r - exp(-1)) < 3 / sqrt(P)
end

@testset "stepper/jitter_overflow" begin
    # an overflowing rate throws instead of giving NaN photons
    function ow(f)
        w = SimWorld(StableRNG(3), cam32(), [jpop(; density=0.0, brightness_jitter=0.3)]; n_sub=1)
        ps = w.pops[1]
        i = _add_emitter!(w, ps, 0.0)
        f(ps, i)
        return w
    end
    @test_throws DomainError SMLMSim.step!(ow((ps, i) -> ps.lj[i] = 800.0), 0.0, 0.01)
    @test_throws DomainError SMLMSim.step!(ow((ps, i) -> (ps.lj[i] = 0.0; ps.γ[i] = Inf)), 0.0, 0.01)
end

@testset "stepper/jitter_zero_alloc" begin
    sig = Population(name=:sig, density=2.0, fluor=two_state(1e4, 20.0, 20.0), lifetime=0.1, budget=5000.0,
                     mobility=[(1.0, 0.1)], brightness_jitter=0.3, jitter_time=0.01, psf=GaussianPSF(0.1))
    oof = Population(name=:oof, layer=:oof, density=0.5, fluor=one_state(5000.0), lifetime=0.1, budget=2e4,
                     z=(0.5, 1.0), brightness_jitter=0.2, psf=GaussianPSF(0.39))
    w = SimWorld(StableRNG(108), cam32(), [sig, oof]; n_sub=4)
    exc = EvanescentExcitation(; depth=0.1, stray=0.05)
    scam = SCMOSCamera(32, 32, 0.1, 1.6; offset=100.0, gain=2.0, qe=1.0)
    nrng = StableRNG(99)
    dst = zeros(32, 32)
    for k in 1:200
        SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01, exc)
    end
    count_allocs(w, exc, nrng, dst, scam, 201)
    total = 0
    for k in 202:301
        total += count_allocs(w, exc, nrng, dst, scam, k)
    end
    @test total == 0
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
    # (a) per emitter height: signal and :oof emitters, excitation, photons and budget draw-down follow I(z)
    sigp = Population(density=20.0, fluor=one_state(1e4), z=(0.0, 0.3), budget=1e9, psf=GaussianPSF(0.1))
    oofp = Population(name=:oof, layer=:oof, density=5.0, fluor=one_state(1e4), z=(0.5, 1.0), budget=1e9,
                      psf=GaussianPSF(0.39))
    w = SimWorld(StableRNG(107), cam32(), [sigp, oofp]; n_sub=1)
    bud0 = Dict((k, ps.id[i]) => ps.budget[i] for (k, ps) in enumerate(w.pops) for i in 1:ps.n)
    SMLMSim.step!(w, 0.0, 0.01, EvanescentExcitation(; depth=0.1, stray=0.05))
    rows = collect(frame_truth(w))
    @test any(r -> w.pops[r.pop].p.layer === :signal && r.z > 0.2, rows)
    @test any(r -> w.pops[r.pop].p.layer === :oof, rows)
    ok = true
    for r in rows
        I = 0.05 + 0.95 * exp(-max(r.z, 0) / 0.1)
        ps = w.pops[r.pop]
        i = findfirst(==(r.id), view(ps.id, 1:ps.n))
        ok &= isapprox(r.excitation, I; rtol=1e-12) && isapprox(r.photons, 1e4 * 0.01 * I; rtol=1e-12) &&
              isapprox(bud0[(Int(r.pop), r.id)] - ps.budget[i], r.photons; rtol=1e-6)
    end
    @test ok
    # (b) a blinking emitter at height z leaves state 1 at 1000 I per second: lit = min(1, 1/(1000 I T))
    fl = GenericFluor(; γ=1e4, q=[-1000.0 1000.0; 1e-9 -1e-9])
    w = SimWorld(StableRNG(109), cam32(), [Population(density=20.0, fluor=fl, z=(0.0, 0.3), psf=GaussianPSF(0.1))]; n_sub=1)
    ps = w.pops[1]
    ps.state[1:ps.n] .= 1
    ps.clock[1:ps.n] .= 1.0
    SMLMSim.step!(w, 0.0, 0.01, EvanescentExcitation(; depth=0.1))
    rows = collect(frame_truth(w))
    ok = true
    for r in rows
        I = exp(-r.z / 0.1)
        lit = min(1.0, 1 / (1000 * I * 0.01))
        ok &= isapprox(r.lit, lit; rtol=1e-9) && isapprox(r.photons, 1e4 * I * 0.01 * lit; rtol=1e-9)
    end
    @test ok
    @test any(r -> r.lit < 1, rows) && any(r -> r.lit == 1, rows)
end
