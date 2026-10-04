using SMLMSim, Test, Statistics, Random, StableRNGs, MicroscopePSFs, Distributions
using SMLMSim.Stepper: _add_emitter!, _reflect, _wrap

# 32x32 pixel camera, 0.1 um pixels: FOV 0.0..3.2 um
cam32() = IdealCamera(1:32, 1:32, 0.1)
one_state(γ) = GenericFluor(; γ, q=zeros(1, 1))
two_state(γ, k_off, k_on) = GenericFluor(; γ, q=[-k_off k_off; k_on -k_on])

# a world with no emitters, ready for _add_emitter!
empty_world(rng, pop; n_sub=1, kw...) = SimWorld(rng, cam32(), [pop]; n_sub, kw...)

function place!(ps, xs, ys)
    for k in eachindex(xs)
        ps.x[k] = xs[k]
        ps.y[k] = ys[k]
    end
end

struct SwitchExc
    ts::Float64
    declared::Bool
end
(e::SwitchExc)(x::Float64, y::Float64, z::Float64, t::Float64) = t < e.ts ? 1.0 : 3.0
SMLMSim.next_switch(e::SwitchExc, t::Float64) = e.declared && t < e.ts ? e.ts : Inf

struct PeriodicExc
    period::Float64
end
(e::PeriodicExc)(x::Float64, y::Float64, z::Float64, t::Float64) = isodd(floor(Int, t / e.period)) ? 2.0 : 1.0
function SMLMSim.next_switch(e::PeriodicExc, t::Float64)
    k = floor(t / e.period) + 1
    return k * e.period > t ? k * e.period : (k + 1) * e.period
end

@testset "stepper/no_distributions" begin
    @test !isdefined(SMLMSim.Stepper, :Distributions)
    dir = joinpath(pkgdir(SMLMSim), "src", "stepper")
    for f in readdir(dir; join=true)
        endswith(f, ".jl") || continue
        for line in eachline(f)
            @test !occursin(r"^\s*(using|import)\s.*Distributions", line)
        end
    end
end

@testset "stepper/energy" begin
    pop = Population(density=1.5, fluor=one_state(1000.0), psf=GaussianPSF(0.13))
    w = SimWorld(StableRNG(1), cam32(), [pop]; n_sub=4, margin=0.0)
    ps = w.pops[1]
    n = ps.n
    @test n > 5
    r = MersenneTwister(2)
    place!(ps, 0.8 .+ 1.6 .* rand(r, n), 0.8 .+ 1.6 .* rand(r, n))
    SMLMSim.step!(w, 0.0, 0.01)
    @test isapprox(sum(w.signal), n * 1000.0 * 0.01; rtol=2e-6)
    @test w.expected === SMLMSim.layers(w).expected
    # generic names stay unexported (SciML and Agents.jl export a step!)
    @test !Base.isexported(SMLMSim, :step!) && !Base.isexported(SMLMSim, :layers)
    @test sum(w.oof) == 0
end

@testset "stepper/bleach_decay" begin
    B = 1000.0 * 0.05
    for I in (1.0, 2.0)
        pop = Population(density=300.0, fluor=one_state(1000.0), budget=B, psf=GaussianPSF(0.13))
        w = SimWorld(StableRNG(3), cam32(), [pop]; n_sub=2)
        N0 = w.pops[1].n
        exc = (x, y, z, t) -> I
        t = 0.0
        for k in 1:10
            SMLMSim.step!(w, t, t + 0.01, exc)
            t += 0.01
            if k in (2, 4, 6, 8, 10)
                p = exp(-1000.0 * I * t / B)
                @test abs(w.pops[1].n / N0 - p) <= 3 * sqrt(p * (1 - p) / N0) + 1e-12
            end
        end
    end
    # a cluster of 10 fluorophores
    pop = Population(density=100.0, fluor=one_state(1000.0), budget=B, multiplicity=10, psf=GaussianPSF(0.13))
    w = SimWorld(StableRNG(4), cam32(), [pop]; n_sub=2)
    N0 = w.pops[1].n
    t = 0.0
    for k in 1:10
        SMLMSim.step!(w, t, t + 0.01)
        t += 0.01
        if k in (2, 5, 10)
            p = exp(-1000.0 * t / B)
            left = sum(w.pops[1].m[1:w.pops[1].n]) / N0
            @test abs(left - 10p) <= 3 * sqrt(10 * p * (1 - p) / N0) + 1e-12
        end
    end
end

@testset "stepper/blink_duty" begin
    k_off, k_on = 100.0, 50.0
    for (I, duty) in ((1.0, k_on / (k_on + k_off)), (2.0, k_on / (k_on + 2k_off)))
        pop = Population(density=400.0, fluor=two_state(1000.0, k_off, k_on), psf=GaussianPSF(0.13))
        w = SimWorld(StableRNG(5), cam32(), [pop]; n_sub=2)
        exc = (x, y, z, t) -> I
        SMLMSim.step!(w, 0.0, 0.3, exc)      # relax to the stationary state at this intensity
        t = 0.3
        fr = Float64[]
        for _ in 1:12
            SMLMSim.step!(w, t, t + 0.1, exc)
            t += 0.1
            ps = w.pops[1]
            push!(fr, count(==(1), ps.state[1:ps.n]) / ps.n)
        end
        se = sqrt(duty * (1 - duty) / w.pops[1].n / length(fr))
        @test abs(mean(fr) - duty) <= 3 * se
    end
end

@testset "stepper/steady_density" begin
    cam = IdealCamera(1:16, 1:16, 0.1)
    pop = Population(density=1.0, lifetime=0.05, fluor=one_state(1000.0), psf=GaussianPSF(0.13))
    w = SimWorld(StableRNG(6), cam, [pop]; n_sub=1)
    A = (w.box[2] - w.box[1]) * (w.box[4] - w.box[3])
    t = 0.0
    counts = Float64[]
    for _ in 1:1500
        SMLMSim.step!(w, t, t + 0.5)
        t += 0.5
        push!(counts, w.pops[1].n)
    end
    @test abs(mean(counts) - A) <= 3 * std(counts) / sqrt(length(counts))
    @test 0.9 <= var(counts) / mean(counts) <= 1.1
end

@testset "stepper/mobility" begin
    cam = IdealCamera(1:64, 1:64, 0.1)
    pop = Population(density=40.0, mobility=[(0.6, 0.5), (0.4, 0.1)], fluor=one_state(1000.0), psf=GaussianPSF(0.13))
    w = SimWorld(StableRNG(7), cam, [pop]; n_sub=4, margin=5.0)
    ps = w.pops[1]
    n = ps.n
    f = count(==(0.5), ps.D[1:n]) / n
    @test abs(f - 0.6) <= 3 * sqrt(0.6 * 0.4 / n)
    x0, y0, D = ps.x[1:n], ps.y[1:n], ps.D[1:n]
    SMLMSim.step!(w, 0.0, 0.1)
    @test ps.n == n
    for Dc in (0.5, 0.1)
        sel = [i for i in 1:n if D[i] == Dc && 2.0 < x0[i] < 9.0 && 2.0 < y0[i] < 9.0]
        msd = mean((ps.x[i] - x0[i])^2 + (ps.y[i] - y0[i])^2 for i in sel)
        @test abs(msd - 4Dc * 0.1) <= 3 * 4Dc * 0.1 / sqrt(length(sel))
    end
end

@testset "stepper/walls" begin
    pop = Population(density=300.0, mobility=[(1.0, 5.0)], fluor=one_state(1000.0), psf=GaussianPSF(0.13))
    for bnd in (:reflecting, :periodic)
        w = SimWorld(StableRNG(8), cam32(), [pop]; n_sub=2, boundary=bnd)
        xmin, xmax, ymin, ymax = w.box
        t = 0.0
        for _ in 1:5
            SMLMSim.step!(w, t, t + 0.5)
            t += 0.5
        end
        ps = w.pops[1]
        n = ps.n
        @test all(xmin .<= ps.x[1:n] .<= xmax) && all(ymin .<= ps.y[1:n] .<= ymax)
        bx = clamp.(ceil.(Int, 10 .* (ps.x[1:n] .- xmin) ./ (xmax - xmin)), 1, 10)
        by = clamp.(ceil.(Int, 10 .* (ps.y[1:n] .- ymin) ./ (ymax - ymin)), 1, 10)
        cnt = zeros(10, 10)
        for k in 1:n
            cnt[bx[k], by[k]] += 1
        end
        χ2 = sum((cnt .- n / 100) .^ 2 ./ (n / 100))
        @test ccdf(Chisq(99), χ2) > 1e-3
    end
    @test _reflect(3.2 + 0.3, 0.0, 3.2) ≈ 2.9
    @test _reflect(-0.3 - 3.2, 0.0, 3.2) ≈ 2.9
    @test _wrap(3.2 + 0.3, 0.0, 3.2) ≈ 0.3
    @test _wrap(-0.3, 0.0, 3.2) ≈ 2.9
end

@testset "stepper/gap_and_order" begin
    B = 1000.0 * 0.05
    pop = Population(density=300.0, fluor=one_state(1000.0), budget=B, psf=GaussianPSF(0.13))
    w = SimWorld(StableRNG(9), cam32(), [pop]; n_sub=4)
    N0 = w.pops[1].n
    SMLMSim.step!(w, 0.0, 0.01)
    @test w.n_gap == 0
    SMLMSim.step!(w, 0.05, 0.06)      # gap of 0.04 s advanced unrecorded
    h = (0.06 - 0.05) / 4
    @test w.n_gap == ceil(Int, (0.05 - 0.01) / h)
    @test w.t == 0.06 && w.frame == 2
    p = exp(-1000.0 * 0.06 / B)
    @test abs(w.pops[1].n / N0 - p) <= 3 * sqrt(p * (1 - p) / N0)
    @test_throws ArgumentError SMLMSim.step!(w, 0.05, 0.07)     # earlier than world time
    @test_throws ArgumentError SMLMSim.step!(w, 0.06, 0.06)     # t_b <= t_a
    @test_throws ArgumentError SMLMSim.step!(w, 0.06, 0.05)
end

@testset "stepper/contiguous_time" begin
    pop = Population(density=0.5, fluor=one_state(1000.0), psf=GaussianPSF(0.13))
    w = SimWorld(StableRNG(10), cam32(), [pop]; n_sub=1)
    T = 0.01
    ok = true
    gaps0 = true
    for k in 1:10_000
        SMLMSim.step!(w, (k - 1) * T, (k - 1) * T + T)
        gaps0 &= w.n_gap == 0
    end
    @test gaps0
    @test w.frame == 10_000
    # with a stretch of 200T the level changes exactly at frames 200j + 1
    bgw = SimWorld(StableRNG(11), cam32(), Population[]; n_sub=1,
                   background=BackgroundModel(level=Uniform(1.0, 100.0), stretch=200T))
    lv = Float64[]
    for k in 1:1000
        SMLMSim.step!(bgw, (k - 1) * T, (k - 1) * T + T)
        push!(lv, bgw.expected[1, 1])
    end
    changes = [k for k in 2:1000 if !isapprox(lv[k], lv[k-1]; rtol=1e-9)]
    @test changes == [201, 401, 601, 801]
end

@testset "stepper/excitation_z" begin
    exc = (x, y, z, t) -> z > 0.25 ? 0.0 : 1.0
    flu = two_state(1000.0, 100.0, 50.0)
    inf = Population(name=:focus, density=1.0, fluor=one_state(1000.0), psf=GaussianPSF(0.13))
    oof = Population(name=:oof, layer=:oof, density=3.0, fluor=flu, budget=50.0, z=(0.5, 1.0),
                     psf=GaussianPSF(0.39))
    w = SimWorld(StableRNG(11), cam32(), [inf, oof]; n_sub=4, margin=0.0)
    ps, po = w.pops
    n = ps.n
    @test n > 3 && po.n > 5
    r = MersenneTwister(12)
    place!(ps, 0.8 .+ 1.6 .* rand(r, n), 0.8 .+ 1.6 .* rand(r, n))
    budget0 = copy(po.budget[1:po.n])
    for k in 1:5
        SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01, exc)
        @test all(iszero, w.oof)
        rows = [r for r in SMLMSim.frame_truth(w) if r.pop == 1]
        @test length(rows) == n
        @test all(r -> r.excitation == 1.0, rows)
        @test isapprox(sum(r.photons for r in rows), sum(w.signal); rtol=2e-6)
    end
    @test po.budget[1:po.n] == budget0
    @test all(isinf, ps.budget[1:n]) && all(==(1), ps.state[1:n])
    @test isapprox(sum(w.signal), n * 1000.0 * 0.01; rtol=2e-6)
end

@testset "stepper/excitation_switch" begin
    pop = Population(density=0.0, fluor=one_state(1000.0), psf=GaussianPSF(0.05))
    ta, tb, ts = 0.0, 0.01, 0.0043
    for declared in (true, false)
        w = SimWorld(StableRNG(13), cam32(), [pop]; n_sub=4, margin=0.0)
        ps = w.pops[1]
        _add_emitter!(w, ps, 0.0)
        place!(ps, [1.6], [1.6])
        SMLMSim.step!(w, ta, tb, SwitchExc(ts, declared))
        te = declared ? ts : 0.005      # undeclared: the switch acts from the next sub-step start
        @test isapprox(sum(w.signal), 1000.0 * (1 * (te - ta) + 3 * (tb - te)); rtol=1e-12)
    end
end

@testset "stepper/birth_energy" begin
    γ, T, life, rate = 5000.0, 0.01, 0.002, 200.0
    pop = Population(density=rate * life, lifetime=life, birth_rate=rate, fluor=one_state(γ),
                     psf=GaussianPSF(0.13))
    w = SimWorld(StableRNG(14), cam32(), [pop]; n_sub=8)
    tot = Float64[]
    for k in 1:2000
        SMLMSim.step!(w, (k - 1) * T, k * T)
        push!(tot, sum(w.signal))
    end
    A_fov = 3.2^2
    @test abs(mean(tot) - rate * life * A_fov * γ * T) <= 3 * std(tot) / sqrt(length(tot))
end

@testset "stepper/validation" begin
    psf = GaussianPSF(0.13)
    f1 = one_state(1000.0)
    @test_throws ArgumentError Population(density=1.0, fluor=f1, psf=psf, layer=:foo)
    @test_throws ArgumentError Population(density=1.0, fluor=f1, psf=psf, mobility=[(0.5, 0.1), (0.4, 0.2)])
    @test_throws ArgumentError Population(density=1.0, fluor=two_state(1000.0, 10.0, 1.0), psf=psf, multiplicity=3)
    q3 = [-1.0 1.0 0.0; 1.0 -2.0 1.0; 0.0 0.0 0.0]       # state 3 absorbing
    @test_throws ArgumentError Population(density=1.0, fluor=GenericFluor(; γ=1000.0, q=q3), psf=psf)
    @test_throws ArgumentError Population(density=1.0, birth_rate=1.0, fluor=f1, psf=psf)
    tbl = StampTable(GaussianPSF(0.13), 0.1, range(-0.3, 0.3, length=5); radius=4)
    @test_throws ArgumentError Population(density=1.0, fluor=f1, psf=tbl, z=(0.5, 1.0))
    @test Population(density=1.0, fluor=f1, psf=tbl, z=(-0.2, 0.2)) isa Population
    @test tbl.pixel_size == 0.1
    p32 = Population(density=1.0, fluor=f1, psf=GaussianPSF(0.13f0))
    @test p32.psf.σ == Float64(0.13f0)
    oofp = Population(layer=:oof, density=1.0, fluor=f1, psf=tbl, z=(-0.2, 0.2))
    @test_throws ArgumentError SimWorld(StableRNG(1), IdealCamera(1:32, 1:32, 0.2), [oofp]; n_sub=1)
    @test SimWorld(StableRNG(1), IdealCamera(1:32, 1:32, 0.1), [oofp]; n_sub=1) isa SimWorld
    good = Population(density=1.0, fluor=f1, psf=psf)
    @test SimWorld(StableRNG(1), cam32(), [good]; n_sub=1, background=BackgroundModel(level=1.0)) isa SimWorld
    @test_throws ArgumentError SimWorld(StableRNG(1), cam32(), [good]; n_sub=1, background=BackgroundModel(level=1.0, stretch=0.0))
    @test_throws ArgumentError SimWorld(StableRNG(1), cam32(), [good]; n_sub=1, boundary=:sticky)
    @test_throws UndefKeywordError SimWorld(StableRNG(1), cam32(), [good])
    @test_throws ArgumentError SimWorld(StableRNG(1), IdealCamera([0.0, 0.1, 0.25, 0.3], [0.0, 0.1, 0.2, 0.3]), [good]; n_sub=1)
    # the states of newborns follow the stationary distribution
    q = [-3.0 2.0 1.0; 1.0 -4.0 3.0; 5.0 1.0 -6.0]
    pop = Population(density=0.0, fluor=GenericFluor(; γ=1000.0, q=q), psf=psf)
    w = SimWorld(StableRNG(15), cam32(), [pop]; n_sub=1)
    ps = w.pops[1]
    for _ in 1:10_000
        _add_emitter!(w, ps, 0.0)
    end
    @test ps.n == 10_000
    π0 = ps.π0
    @test isapprox(π0' * q, zeros(1, 3); atol=1e-12) && sum(π0) ≈ 1
    obs = [count(==(UInt8(s)), ps.state[1:ps.n]) for s in 1:3]
    χ2 = sum((obs .- 10_000 .* π0) .^ 2 ./ (10_000 .* π0))
    @test ccdf(Chisq(2), χ2) > 1e-3
end

const GOLDEN_SUM = 54694.319025331104

function rich_world(seed)
    flu = two_state(3000.0, 200.0, 20.0)
    tbl = StampTable(GaussianPSF(0.39), 0.1, range(-0.3, 1.0, length=8); radius=12)
    sig = Population(name=:sig, density=1.0, lifetime=0.1, fluor=flu, budget=100.0,
                     mobility=[(0.5, 0.3), (0.5, 0.05)], psf=GaussianPSF(0.13))
    oof = Population(name=:oof, layer=:oof, density=1.0, lifetime=0.05, fluor=flu, budget=200.0,
                     z=(0.5, 1.0), psf=tbl)
    oof2 = Population(name=:oof2, layer=:oof, density=1.0, lifetime=0.05, fluor=flu, budget=200.0,
                      z=(-0.2, 0.2), psf=GaussianPSF(0.39))
    bg = BackgroundModel(level=5000.0, jitter=0.05, contrast=0.3, feature_size=0.3,
                         correlation_time=0.05, illumination_width=2.0, stretch=0.1)
    return SimWorld(StableRNG(seed), cam32(), [sig, oof, oof2]; n_sub=4, background=bg, merge_radius=0.25)
end

@testset "stepper/determinism" begin
    a, b = rich_world(21), rich_world(21)
    for k in 1:50
        ea = copy(SMLMSim.step!(a, (k - 1) * 0.01, k * 0.01))
        eb = SMLMSim.step!(b, (k - 1) * 0.01, k * 0.01)
        @test ea == eb
        @test collect(SMLMSim.frame_truth(a)) == collect(SMLMSim.frame_truth(b))
    end
    @test a.pops[1].x[1:a.pops[1].n] == b.pops[1].x[1:b.pops[1].n]
    g = rich_world(22)
    for k in 1:20
        SMLMSim.step!(g, (k - 1) * 0.01, k * 0.01)
    end
    @test isapprox(sum(g.expected), GOLDEN_SUM; rtol=1e-12)
end

const_excitation(x, y, z, t) = 2.0

function count_allocs(w, exc::E, nrng, dst, cam, k) where {E}
    return @allocated begin
        SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01, exc)
        SMLMSim.scmos_noise!(nrng, copyto!(dst, w.expected), cam)
    end
end

@testset "stepper/zero_alloc" begin
    w = rich_world(23)
    exc = PeriodicExc(0.0173)
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
    # a closure and a named function specialise like a functor
    for exc2 in (let I = 2.0; (x, y, z, t) -> I end, const_excitation)
        w2 = rich_world(23)
        for k in 1:200
            SMLMSim.step!(w2, (k - 1) * 0.01, k * 0.01, exc2)
        end
        count_allocs(w2, exc2, nrng, dst, scam, 201)
        total = 0
        for k in 202:301
            total += count_allocs(w2, exc2, nrng, dst, scam, k)
        end
        @test total == 0
    end
end

@testset "stepper/zero_alloc_level_draw" begin
    # a distribution level redrawn every 5 frames allocates nothing either
    w = SimWorld(StableRNG(24), cam32(), Population[]; n_sub=1,
                 background=BackgroundModel(level=Uniform(1.0, 100.0), stretch=0.05, contrast=0.3))
    for k in 1:20
        SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01)
    end
    total = 0
    for k in 21:120
        total += @allocated SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01)
    end
    @test total == 0
end

@testset "stepper/background_gap_draws" begin
    # every stretch boundary crossed draws a level, including those inside a gap
    mk() = SimWorld(StableRNG(123), cam32(), Population[]; n_sub=1,
                    background=BackgroundModel(level=Uniform(1.0, 100.0), stretch=1.0))
    a = mk(); for t in (1.0, 2.0, 3.0); SMLMSim.step!(a, t, t + 0.01); end
    b = mk(); SMLMSim.step!(b, 3.0, 3.01)
    r = StableRNG(123); L = [rand(r, Uniform(1.0, 100.0)) for _ in 1:4]
    @test a.bg.level == b.bg.level == L[4]          # 5.230326015385788
    @test rand(a.rng) == rand(b.rng)
end

@testset "stepper/nonfinite_inputs" begin
    psf = GaussianPSF(0.13)
    mkpop(; kw...) = Population(; density=1.0, fluor=one_state(1000.0), psf, kw...)
    mkworld(; kw...) = SimWorld(StableRNG(1), cam32(), [mkpop()]; n_sub=1, kw...)
    for exc in ((x, y, z, t) -> NaN, (x, y, z, t) -> Inf, (x, y, z, t) -> -1.0)
        w = mkworld()
        @test w.pops[1].n > 0
        @test_throws DomainError SMLMSim.step!(w, 0.0, 0.01, exc)
    end
    @test_throws ArgumentError mkpop(density=Inf)
    @test_throws ArgumentError mkpop(birth_rate=Inf, lifetime=1.0)
    @test_throws ArgumentError mkpop(brightness_sigma=Inf)
    @test_throws ArgumentError mkpop(mobility=[(1.0, Inf)])
    @test_throws ArgumentError mkpop(z=(0.0, Inf))
    @test_throws ArgumentError mkpop(fluor=one_state(Inf))
    @test_throws ArgumentError mkpop(fluor=one_state(NaN))
    @test_throws ArgumentError mkpop(fluor=GenericFluor(; γ=1000.0, q=[-Inf Inf; 1.0 -1.0]))
    @test_throws ArgumentError mkworld(margin=Inf)
    @test_throws ArgumentError mkworld(t0=NaN)
    @test_throws ArgumentError SimWorld(StableRNG(1), cam32(), [mkpop(density=1e308)]; n_sub=1)
    @test_throws ArgumentError mkworld(background=BackgroundModel(level=Inf))
    @test_throws ArgumentError mkworld(background=BackgroundModel(jitter=Inf))
    @test_throws ArgumentError mkworld(background=BackgroundModel(contrast=Inf))
    @test_throws DomainError mkworld(background=BackgroundModel(level=Uniform(-2.0, -1.0)))
    @test_throws ArgumentError SMLMSim.step!(mkworld(), 0.0, Inf)
    @test_throws ArgumentError SMLMSim.step!(mkworld(), NaN, 0.01)
    @test mkworld() isa SimWorld
end

@testset "stepper/finite_inputs_converted" begin
    psf = GaussianPSF(0.13)
    mkpop(; kw...) = Population(; density=1.0, fluor=one_state(1000.0), psf, kw...)
    mkworld(; kw...) = SimWorld(StableRNG(1), cam32(), [mkpop()]; n_sub=1, kw...)
    # Codex's input: every row passes the sum test with an infinite off-diagonal
    @test_throws ArgumentError mkpop(fluor=GenericFluor(; γ=1000.0, q=[-1.0 Inf; 1.0 -1.0]), multiplicity=0)
    @test_throws ArgumentError mkpop(fluor=GenericFluor(; γ=1000.0, q=[-1.0 1.0; NaN -1.0]), multiplicity=0)
    # a number that overflows Float64 is checked after the conversion
    @test_throws ArgumentError mkworld(background=BackgroundModel(level=big"1e400"))
    @test_throws ArgumentError mkworld(t0=big"1e400")
    @test_throws ArgumentError mkworld(margin=big"1e400")
    w = mkworld(t0=big"0.5")
    @test w.t === 0.5
end

# Every number a stepper input type takes is checked as the Float64 it is stored as. The walk below takes each
# numeric input of each type, hands it as a BigFloat beyond or below the Float64 range, and compares the outcome
# with the Float64 conversion of the same number: both throw, or both construct with isequal stored values. The
# field lists are reflected, so a new numeric field fails here until a case covers it.
@testset "stepper/inputs_as_stored" begin
    bigs = (big"1e400", big"1e-400", -big"1e400", -big"1e-400")
    snap(x::Union{Number,Symbol,Bool}) = x
    snap(x::AbstractArray) = map(snap, x)
    snap(x::Tuple) = map(snap, x)
    snap(x) = map(f -> snap(getfield(x, f)), fieldnames(typeof(x)))      # structs: field values, not types
    # any other exception (a MethodError, an InexactError) is an outcome of its own and fails the comparison
    attempt(f) = try (:ok, f()) catch e; (e isa Union{ArgumentError,DomainError} ? :throw : :error, typeof(e)) end

    # a case varies one input: `mk(v)` builds with v and returns what the world stores; `fin` is true when the docs
    # require a finite value, `pos` when they require > 0
    cases = NamedTuple{(:T, :field, :label, :mk, :fin, :pos),Tuple{Symbol,Union{Symbol,Nothing},String,Function,Bool,Bool}}[]
    add!(T, field, label, mk; fin=false, pos=false) = push!(cases, (; T, field, label, mk, fin, pos))

    # ---- Population: every keyword that is a Float64, or a tuple or vector of them, then fluor and psf
    g0, q0 = 1000.0, [-5.0 5.0; 2.0 -2.0]
    pkw = (density=1.0, lifetime=5.0, birth_rate=0.2, mobility=[(1.0, 0.1), (0.0, 0.2)], fluor=GenericFluor(; γ=g0, q=q0),
           brightness_sigma=0.1, budget=1e5, multiplicity=1, z=(0.0, 0.1), psf=GaussianPSF(0.13),
           brightness_jitter=0.1, jitter_time=0.01)
    function pop_snap(; kw...)
        p = Population(; merge(pkw, kw)...)
        w = SimWorld(StableRNG(1), cam32(), [p]; n_sub=1)
        ps = w.pops[1]
        return (map(f -> snap(getfield(p, f)), filter(!=(:fluor), collect(fieldnames(Population)))),
                ps.γ0, ps.q, ps.exitrate, ps.sigma_px, ps.n, w.box)
    end
    for (k, fin, pos) in ((:density, true, false), (:lifetime, false, true), (:birth_rate, true, false),
                          (:brightness_sigma, true, false), (:budget, false, true),
                          (:brightness_jitter, true, false), (:jitter_time, false, true))
        add!(:Population, k, string(k), v -> pop_snap(; NamedTuple{(k,)}((v,))...); fin, pos)
    end
    add!(:Population, :mobility, "mobility fraction 1", v -> pop_snap(; mobility=[(v, 0.1), (0.0, 0.2)]))
    add!(:Population, :mobility, "mobility D 1", v -> pop_snap(; mobility=[(1.0, v), (0.0, 0.2)]); fin=true)
    add!(:Population, :mobility, "mobility fraction 2", v -> pop_snap(; mobility=[(1.0, 0.1), (v, 0.2)]))
    add!(:Population, :mobility, "mobility D 2", v -> pop_snap(; mobility=[(1.0, 0.1), (0.0, v)]); fin=true)
    add!(:Population, :z, "z 1", v -> pop_snap(; z=(v, 0.1)); fin=true)
    add!(:Population, :z, "z 2", v -> pop_snap(; z=(0.0, v)); fin=true)
    add!(:Population, nothing, "fluor.γ", v -> pop_snap(; fluor=GenericFluor(; γ=v, q=q0)); fin=true)
    add!(:Population, nothing, "psf σ", v -> pop_snap(; psf=GaussianPSF(v)); fin=true, pos=true)
    for (i, j) in ((1, 2), (2, 1))              # one off-diagonal entry, its diagonal set so the row sums to 0 in BigFloat
        function qmk(v)
            q = big.(q0)
            q[i, j] = v
            q[i, i] = -v
            return pop_snap(; fluor=GenericFluor(; γ=g0, q))
        end
        add!(:Population, nothing, "fluor.q[$i,$j]", qmk; fin=true)
    end
    # StampTable: pixel_size and the two ends of zs
    stamp_snap(px, zs) = (t = StampTable(GaussianPSF(0.13), px, zs; radius=3, oversample=2);
                          (t.stamps, collect(t.zs), t.radius, t.oversample, t.zinterp, t.pixel_size))
    add!(:Population, nothing, "StampTable pixel_size", v -> stamp_snap(v, range(big"-0.3", big"0.3"; length=5)); fin=true, pos=true)
    add!(:Population, nothing, "StampTable zs start", v -> stamp_snap(big"0.1", range(v, big"0.3"; length=5)); fin=true)
    add!(:Population, nothing, "StampTable zs stop", v -> stamp_snap(big"0.1", range(big"-0.3", v; length=5)); fin=true)

    # ---- DimerKinetics: every field; D_dimer is a number here
    dkw = (k_on=10.0, r_react=0.03, k_off=0.2, D_rot=1.0, d_dimer=0.01, D_dimer=0.05)
    dim(; kw...) = snap(DimerKinetics(; merge(dkw, kw)...))
    for (k, fin, pos) in ((:k_on, false, true), (:r_react, true, true), (:k_off, true, false), (:D_rot, true, false),
                          (:d_dimer, true, false), (:D_dimer, true, false))
        add!(:DimerKinetics, k, string(k), v -> dim(; NamedTuple{(k,)}((v,))...); fin, pos)
    end

    # ---- BackgroundModel: every field, through SimWorld, where it is checked
    bkw = (level=5.0, stretch=1.0, jitter=0.1, feature_size=0.8, contrast=0.3, correlation_time=2.0, illumination_width=3.0)
    function bg_snap(; kw...)
        w = SimWorld(StableRNG(2), cam32(), [Population(; pkw..., density=0.0, birth_rate=0.0)]; n_sub=1,
                     background=BackgroundModel(; merge(bkw, kw)...))
        b = w.bg
        return (b.level, b.stretch, b.jitter, b.contrast, b.tau, b.P, b.g, b.iy, b.wx)
    end
    for (k, fin, pos) in ((:level, true, false), (:stretch, false, true), (:jitter, true, false),
                          (:feature_size, false, true), (:contrast, true, false), (:correlation_time, false, true),
                          (:illumination_width, false, true))
        add!(:BackgroundModel, k, string(k), v -> bg_snap(; NamedTuple{(k,)}((v,))...); fin, pos)
    end

    # ---- Spot and SpotExcitation (tan_cos and tan_sin are derived from tilt)
    skw = (x=1.0, y=1.5, σ=0.3, gain=2.0, z_R=Inf, t_on=0.1, t_off=0.2)
    for (k, fin, pos) in ((:x, true, false), (:y, true, false), (:σ, true, true), (:gain, true, false),
                          (:z_R, false, true), (:t_on, false, false), (:t_off, false, false))
        add!(:Spot, k, string(k), v -> snap(Spot(; merge(skw, NamedTuple{(k,)}((v,)))...)); fin, pos)
    end
    xkw = (base=1.0, f_evan=0.3, d_evan=0.1, λ=0.642, n=1.33, tilt=(1.0, 0.5))
    ex(; kw...) = snap(SpotExcitation(; spots=[Spot(; skw...)], merge(xkw, kw)...))
    for (k, fin, pos) in ((:base, true, false), (:f_evan, true, false), (:d_evan, false, true),
                          (:λ, true, true), (:n, true, true))
        add!(:SpotExcitation, k, string(k), v -> ex(; NamedTuple{(k,)}((v,))...); fin, pos)
    end
    add!(:SpotExcitation, :tilt, "tilt 1", v -> ex(; tilt=(v, 0.5)); fin=true)
    add!(:SpotExcitation, :tilt, "tilt 2", v -> ex(; tilt=(1.0, v)); fin=true)

    # ---- EvanescentExcitation (Core type, shared with the diffusion path)
    add!(:EvanescentExcitation, :depth, "depth", v -> snap(EvanescentExcitation(; depth=v, stray=0.1)); pos=true)
    add!(:EvanescentExcitation, :stray, "stray", v -> snap(EvanescentExcitation(; depth=0.1, stray=v)); fin=true)

    # ---- SimWorld keywords and step!'s times, listed explicitly (they are keywords, not fields)
    wpop = Population(; pkw...)
    function world_snap(; kw...)
        w = SimWorld(StableRNG(3), cam32(), [wpop]; n_sub=2, kw...)
        return (w.box, w.t, w.t_a, w.t_b, w.merge_radius, w.pops[1].n)
    end
    add!(:SimWorld, :margin, "margin", v -> world_snap(; margin=v); fin=true)
    add!(:SimWorld, :t0, "t0", v -> world_snap(; t0=v); fin=true)
    add!(:SimWorld, :merge_radius, "merge_radius", v -> world_snap(; merge_radius=v))
    function step_snap(ta, tb)
        w = SimWorld(StableRNG(3), cam32(), [wpop]; n_sub=2)
        img = copy(SMLMSim.step!(w, ta, tb))
        return (img, snap(collect(SMLMSim.frame_truth(w))), w.t, w.frame)
    end
    add!(:step!, :t_a, "t_a", v -> step_snap(v, 0.01); fin=true)
    add!(:step!, :t_b, "t_b", v -> step_snap(0.0, v); fin=true)

    # ---- coverage: the numeric fields of each struct, less the derived ones, are exactly the fields varied above
    hasnum(ft) = ft === Float64 || ft isa TypeVar || (ft isa Union ? any(hasnum, Base.uniontypes(ft)) :
                 ft <: Tuple ? any(hasnum, fieldtypes(ft)) : ft <: AbstractVector && hasnum(eltype(ft)))
    numeric(T, derived) = Set(f for (f, ft) in zip(fieldnames(T), fieldtypes(T)) if hasnum(ft) && !(f in derived))
    covered(T) = Set(c.field for c in cases if c.T === T && c.field !== nothing)
    # derived: tan_cos and tan_sin follow from tilt. multiplicity is an Int (not converted), name, layer, binds, spots,
    # fluor and psf are not Float64 inputs; fluor's and psf's numbers have their own cases above.
    @test numeric(Population, ()) == covered(:Population)
    @test numeric(DimerKinetics, ()) == covered(:DimerKinetics)
    @test numeric(BackgroundModel{Float64}, ()) == covered(:BackgroundModel)
    @test numeric(Spot, ()) == covered(:Spot)
    @test numeric(SpotExcitation, (:tan_cos, :tan_sin)) == covered(:SpotExcitation)
    @test numeric(EvanescentExcitation, ()) == covered(:EvanescentExcitation)
    @test covered(:SimWorld) == Set([:margin, :t0, :merge_radius])
    @test covered(:step!) == Set([:t_a, :t_b])
    @test length(cases) == 54

    outcome(r) = r[1] === :ok ? "constructs" : "$(r[1]) $(r[2])"
    nboth = 0                           # inputs both routes construct with (the stored values then compared)
    bad = String[]                      # the failing inputs, named, so the test's failure message lists them
    for c in cases, v in bigs
        a, b = attempt(() -> c.mk(v)), attempt(() -> c.mk(Float64(v)))
        what = "$(c.T) $(c.label), v = $(Float64(v)): BigFloat $(outcome(a)), Float64 $(outcome(b))"
        if a[1] !== b[1] || a[1] === :error
            push!(bad, what)
        elseif a[1] === :ok ? !isequal(a[2], b[2]) : a[2] !== b[2]
            push!(bad, what * " with different stored values")
        end
        nboth += a[1] === b[1] === :ok
    end
    for c in cases
        c.fin && attempt(() -> c.mk(big"1e400"))[1] !== :throw && push!(bad, "$(c.T) $(c.label) accepts 1e400")
        c.pos && attempt(() -> c.mk(big"1e-400"))[1] !== :throw && push!(bad, "$(c.T) $(c.label) accepts 1e-400")
    end
    @test bad == String[]
    @test nboth >= 30

    # ---- Codex's round-3 inputs, by name
    @test_throws ArgumentError Population(density=1.0, multiplicity=0, fluor=GenericFluor(; γ=big"1e400", q=zeros(1, 1)),
                                          psf=GaussianPSF(0.13))
    @test_throws ArgumentError Population(density=0.0, birth_rate=1.0, budget=1.0, lifetime=Inf,
                                          fluor=GenericFluor(; γ=big"1e-400", q=zeros(1, 1)), psf=GaussianPSF(0.13))
    M = BigFloat(floatmax(Float64))
    qM = [-M M/3 M/3 M/3; 1 -1 0 0; 1 0 -1 0; 1 0 0 -1]
    @test_throws ArgumentError Population(density=1.0, multiplicity=0, fluor=GenericFluor(; γ=1000.0, q=qM), psf=GaussianPSF(0.13))
    @test_throws ArgumentError Population(density=0.0, fluor=GenericFluor(; γ=1000.0, q=zeros(1, 1)), psf=GaussianPSF(big"1e-400"))
end

@testset "stepper/ctmc_state_limit" begin
    function cyc(n)
        q = zeros(n, n)
        for i in 1:n
            q[i, mod1(i + 1, n)] = 1.0
            q[i, i] = -1.0
        end
        return q
    end
    mk(n) = Population(density=1.0, fluor=GenericFluor(; γ=1000.0, q=cyc(n)), psf=GaussianPSF(0.13))
    @test_throws ArgumentError mk(256)
    w = SimWorld(StableRNG(5), cam32(), [mk(255)]; n_sub=1)
    @test SMLMSim.step!(w, 0.0, 0.01) isa Matrix{Float64}
end

@testset "stepper/switch_budget_clock" begin
    # budgets and CTMC clocks carry across declared switches and across gaps (coverage)
    function sw_world(fluor; budget=Inf)
        pop = Population(density=0.0, fluor=fluor, budget=budget, psf=GaussianPSF(0.13))
        w = SimWorld(StableRNG(1), cam32(), [pop]; n_sub=4, margin=0.0)
        ps = w.pops[1]
        _add_emitter!(w, ps, 0.0)
        place!(ps, [1.6], [1.6])
        return w, ps
    end
    blink = GenericFluor(; γ=1000.0, q=[-100.0 100.0; 1e-12 -1e-12])
    # (a) the budget runs out after the switch: 4 photons at 1000/s, then 6 at 3000/s
    w, ps = sw_world(one_state(1000.0); budget=1e6)
    ps.budget[1] = 10.0
    SMLMSim.step!(w, 0.0, 0.01, SwitchExc(0.004, true))
    r = only(SMLMSim.frame_truth(w))
    @test r.t_bleach ≈ 0.006 rtol = 1e-9
    @test r.photons ≈ 10 rtol = 1e-9
    @test r.m == 0
    # (b) the clock runs out after the switch
    w, ps = sw_world(blink)
    ps.state[1] = 1; ps.clock[1] = 1.0
    SMLMSim.step!(w, 0.0, 0.01, SwitchExc(0.004, true))
    r = only(SMLMSim.frame_truth(w))
    @test r.lit ≈ 0.6 rtol = 1e-9
    @test r.photons ≈ 10 rtol = 1e-9
    # (c) a switch inside the gap, the clock carried through
    w, ps = sw_world(blink)
    ps.state[1] = 1; ps.clock[1] = 3.5
    SMLMSim.step!(w, 0.0, 0.01, SwitchExc(0.015, true))
    r = only(SMLMSim.frame_truth(w))
    @test r.photons ≈ 10 rtol = 1e-9
    @test r.lit ≈ 1 rtol = 1e-9
    SMLMSim.step!(w, 0.02, 0.03, SwitchExc(0.015, true))
    r = only(SMLMSim.frame_truth(w))
    @test r.photons ≈ 5 rtol = 1e-9
    @test r.lit ≈ 1 / 6 rtol = 1e-9
    # (d) a switch inside the gap, the budget carried through
    w, ps = sw_world(one_state(1000.0); budget=1e6)
    ps.budget[1] = 32.0
    SMLMSim.step!(w, 0.0, 0.01, SwitchExc(0.015, true))
    @test only(SMLMSim.frame_truth(w)).photons ≈ 10 rtol = 1e-9
    SMLMSim.step!(w, 0.02, 0.03, SwitchExc(0.015, true))
    r = only(SMLMSim.frame_truth(w))
    @test r.photons ≈ 2 rtol = 1e-9
    @test r.t_bleach ≈ 0.02 + 1 / 1500 rtol = 1e-9
    @test r.m == 0
end

@testset "stepper/mixed_mobility" begin
    pop = Population(density=40.0, mobility=[(0.5, 0.0), (0.5, 0.3)], fluor=one_state(1000.0),
                     psf=GaussianPSF(0.13))
    for boundary in (:reflecting, :periodic)
        w = SimWorld(StableRNG(7), IdealCamera(1:64, 1:64, 0.1), [pop]; n_sub=4, margin=5.0, boundary)
        ps = w.pops[1]
        n = ps.n
        start = Dict(ps.id[i] => (ps.x[i], ps.y[i]) for i in 1:n)
        for k in 1:10
            SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01)
        end
        @test ps.n == n
        still = [(ps.x[i], ps.y[i]) == start[ps.id[i]] for i in 1:n if ps.D[i] == 0.0]
        moved = [(ps.x[i], ps.y[i]) != start[ps.id[i]] for i in 1:n if ps.D[i] == 0.3]
        @test all(still)
        @test count(moved) > 0.95 * length(moved)
        @test abs(length(still) - n / 2) < 3 * sqrt(n / 4)
    end
end
