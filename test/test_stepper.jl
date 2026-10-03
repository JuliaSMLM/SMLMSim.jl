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
    return SimWorld(StableRNG(seed), cam32(), [sig, oof, oof2]; n_sub=4, background=bg)
end

@testset "stepper/determinism" begin
    a, b = rich_world(21), rich_world(21)
    for k in 1:50
        ea = copy(SMLMSim.step!(a, (k - 1) * 0.01, k * 0.01))
        eb = SMLMSim.step!(b, (k - 1) * 0.01, k * 0.01)
        @test ea == eb
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
