using SMLMSim, Test, Statistics, Random, StableRNGs, MicroscopePSFs
using SMLMSim: frame_truth, DimerKinetics, DiffusionSMLMConfig, FrameTruth
using SMLMSim.Stepper: _put_pair!, _add_emitter!, _substep!, _begin_truth!, _finish_truth!, _advance_acc!, _write_row!
using Random: randexp
using Distributions: Uniform
using SMLMSim.InteractionDiffusion: _disp, _norm

# the design section 7a configuration, defined once in the benchmark (included, it runs nothing)
include(joinpath(@__DIR__, "..", "dev", "benchmark_stepper.jl"))

cam32() = IdealCamera(1:32, 1:32, 0.1)
one_state(γ) = GenericFluor(; γ, q=zeros(1, 1))
two_state(γ, k_off, k_on) = GenericFluor(; γ, q=[-k_off k_off; k_on -k_on])

function place!(ps, xs, ys)
    for k in eachindex(xs)
        ps.x[k] = xs[k]
        ps.y[k] = ys[k]
    end
end

# largest distance between the empirical CDF of ts and the exponential CDF of the given mean
function ks_exp(ts, mean)
    s = sort(ts)
    n = length(s)
    return maximum(max(abs(k / n - (1 - exp(-s[k] / mean))), abs((k - 1) / n - (1 - exp(-s[k] / mean))))
                   for k in 1:n)
end

@testset "closedloop/truth_energy" begin
    pop = Population(density=3.0, fluor=two_state(1e4, 300.0, 200.0), budget=3000.0, psf=GaussianPSF(0.05))
    w = SimWorld(StableRNG(1), cam32(), [pop]; n_sub=4, margin=0.0)
    ps = w.pops[1]
    r = MersenneTwister(2)
    place!(ps, 0.4 .+ 2.4 .* rand(r, ps.n), 0.4 .+ 2.4 .* rand(r, ps.n))
    T = 0.02
    for k in 1:40
        SMLMSim.step!(w, (k - 1) * T, k * T)
        rows = frame_truth(w)
        @test isapprox(sum(r.photons for r in rows if w.pops[r.pop].p.layer === :signal), sum(w.signal); rtol=2e-6)
    end
end

@testset "closedloop/truth_consistency" begin
    # the photon identity holds here because budget = Inf (no bleach) and there is no jitter
    pop = Population(density=3.0, fluor=two_state(1e4, 300.0, 200.0), psf=GaussianPSF(0.05))
    w = SimWorld(StableRNG(3), cam32(), [pop]; n_sub=4, margin=0.0)
    I0 = 1.5
    exc = (x, y, z, t) -> I0
    T = 0.02
    for k in 1:20
        SMLMSim.step!(w, (k - 1) * T, k * T, exc)
        for r in frame_truth(w)
            @test isapprox(r.photons, r.m * 1e4 * I0 * r.lit * T; rtol=1e-12)
            @test isapprox(r.excitation, I0; rtol=1e-12)
            @test 0 <= r.lit <= 1
        end
    end
end

@testset "closedloop/bleach_time" begin
    B, γ = 2000.0, 1e5
    for I0 in (1.0, 2.0)
        pop = Population(density=500.0, fluor=one_state(γ), budget=B, psf=GaussianPSF(0.05))
        w = SimWorld(StableRNG(4), cam32(), [pop]; n_sub=1, margin=0.0)
        SMLMSim.step!(w, 0.0, 1.0, (x, y, z, t) -> I0)
        tb = [r.t_bleach for r in frame_truth(w) if !isnan(r.t_bleach)]
        @test length(tb) > 4000
        @test ks_exp(tb, B / (γ * I0)) * sqrt(length(tb)) < 1.63   # p > 0.01
    end
end

@testset "closedloop/truth_rows" begin
    pop = Population(density=200.0, lifetime=0.02, budget=2000.0, fluor=one_state(1e5), psf=GaussianPSF(0.05))
    w = SimWorld(StableRNG(5), cam32(), [pop]; n_sub=4, margin=0.0)
    T = 0.01
    ended = Set{Int}()
    for k in 1:100
        ta, tb = (k - 1) * T, k * T
        npre = k == 1 ? 0 : sum(ps.n for ps in w.pops)   # initial emitters are born at t0 = t_a
        SMLMSim.step!(w, ta, tb)
        rows = frame_truth(w)
        ids = [r.id for r in rows]
        @test allunique(ids)
        @test length(rows) == npre + count(r -> !isnan(r.t_birth), rows)
        @test isdisjoint(ids, ended)
        for r in rows
            (isnan(r.t_depart) && isnan(r.t_bleach)) || push!(ended, r.id)
        end
        @test count(r -> !all(t -> isnan(t) || (ta <= t < tb), (r.t_birth, r.t_depart, r.t_bleach)), rows) == 0
        @test count(r -> !(0 <= r.lit <= 1), rows) == 0
        @test count(r -> r.frame != k, rows) == 0
    end
    @test !isempty(ended)
end

@testset "closedloop/overlap" begin
    flag(w) = Dict(r.id => r.overlap for r in frame_truth(w))
    pop = Population(density=0.0, fluor=one_state(1000.0), psf=GaussianPSF(0.05))
    w = SimWorld(StableRNG(6), cam32(), [pop]; n_sub=1, margin=0.0, merge_radius=0.25)
    ps = w.pops[1]
    for _ in 1:3
        _add_emitter!(w, ps, 0.0)
    end
    place!(ps, [1.0, 1.2, 1.6], [1.0, 1.0, 1.0])
    SMLMSim.step!(w, 0.0, 0.01)
    f = flag(w)
    @test f[ps.id[1]] && f[ps.id[2]] && !f[ps.id[3]]
    # the flag is recomputed every exposure: move the second emitter 0.6 μm from the first and none is flagged
    ps.x[2], ps.y[2] = 1.0, 1.6
    SMLMSim.step!(w, 0.01, 0.02)
    @test !any(values(flag(w)))
    w0 = SimWorld(StableRNG(6), cam32(), [pop]; n_sub=1, margin=0.0)
    for _ in 1:2
        _add_emitter!(w0, w0.pops[1], 0.0)
    end
    place!(w0.pops[1], [1.0, 1.2], [1.0, 1.0])
    SMLMSim.step!(w0, 0.0, 0.01)
    @test !any(values(flag(w0)))
    wp = SimWorld(StableRNG(6), cam32(), [pop]; n_sub=1, margin=0.0, merge_radius=0.25, boundary=:periodic)
    for _ in 1:2
        _add_emitter!(wp, wp.pops[1], 0.0)
    end
    place!(wp.pops[1], [0.05, 3.15], [1.0, 1.0])
    SMLMSim.step!(wp, 0.0, 0.01)
    @test all(values(flag(wp)))
    oof = Population(name=:oof, layer=:oof, density=0.0, fluor=one_state(1000.0), psf=GaussianPSF(0.2))
    wo = SimWorld(StableRNG(6), cam32(), [pop, oof]; n_sub=1, margin=0.0, merge_radius=0.25)
    for ps in wo.pops
        _add_emitter!(wo, ps, 0.0)
    end
    place!(wo.pops[1], [1.0], [1.0])
    place!(wo.pops[2], [1.1], [1.0])
    SMLMSim.step!(wo, 0.0, 0.01)
    @test !any(values(flag(wo)))
    @test_throws ArgumentError SimWorld(StableRNG(6), cam32(), [pop]; n_sub=1, merge_radius=-1)
end

# ---- PR-3 lane B: cross-population dimers ----

struct PeriodicExc
    period::Float64
end
(e::PeriodicExc)(x::Float64, y::Float64, z::Float64, t::Float64) = isodd(floor(Int, t / e.period)) ? 2.0 : 1.0
function SMLMSim.next_switch(e::PeriodicExc, t::Float64)
    k = floor(t / e.period) + 1
    return k * e.period > t ? k * e.period : (k + 1) * e.period
end

byid(rows) = Dict(r.id => r for r in rows)

# a world of two immobile emitters of one population (or one per population) 0.004 um apart
function pair_world(pops, dimers; seed=1, boundary=:reflecting, xs=[1.6, 1.604], ys=[1.6, 1.6], n_sub=4)
    w = SimWorld(StableRNG(seed), cam32(), pops; n_sub, margin=0.0, dimers, boundary)
    if length(pops) == 1
        _add_emitter!(w, w.pops[1], 0.0)
        _add_emitter!(w, w.pops[1], 0.0)
        place!(w.pops[1], xs, ys)
    else
        for (k, ps) in enumerate(w.pops)
            _add_emitter!(w, ps, 0.0)
            place!(ps, xs[k:k], ys[k:k])
        end
    end
    return w
end

@testset "closedloop/dimer_equilibrium" begin
    h, k_on, k_off, r = 1.25e-3, 2000.0, 20.0, 0.01
    pop = Population(density=30.0, fluor=one_state(1000.0), mobility=[(1.0, 1.0)], psf=GaussianPSF(0.1))
    w = SimWorld(StableRNG(11), cam32(), [pop]; n_sub=1, margin=0.0, boundary=:periodic,
                 dimers=DimerKinetics(k_on=k_on, r_react=r, k_off=k_off, D_rot=1.0, d_dimer=0.01))
    A = 3.2^2
    ps = w.pops[1]
    for k in 1:500
        SMLMSim.step!(w, (k - 1) * h, k * h)
    end
    nd = Float64[]
    nf2 = Float64[]
    for k in 501:30500
        SMLMSim.step!(w, (k - 1) * h, k * h)
        d = count(!=(0), view(ps.partner, 1:ps.n)) / 2
        push!(nd, d)
        nf = ps.n - 2d
        push!(nf2, nf * (nf - 1))
    end
    k̃ = (1 - exp(-k_on * h)) / h
    lhs = k_off * mean(nd)
    rhs = k̃ * π * r^2 * mean(nf2) / (2A)
    @test isapprox(lhs, rhs; rtol=0.05)
end

function cross_world(seed, DB; k_on=50.0, binds_B=true, d_dimer=0.05)
    cam = IdealCamera(1:64, 1:64, 0.1)
    A = Population(name=:A, density=20.0, mobility=[(1.0, 4.0)], fluor=one_state(1000.0), psf=GaussianPSF(0.1))
    B = Population(name=:B, density=20.0, mobility=[(1.0, DB)], fluor=one_state(1000.0), psf=GaussianPSF(0.1),
                   binds=binds_B)
    return SimWorld(StableRNG(seed), cam, [A, B]; n_sub=8, margin=0.0, boundary=:periodic,
                    dimers=DimerKinetics(k_on=k_on, r_react=0.01, k_off=20.0, D_rot=1.0, d_dimer=d_dimer))
end

@testset "closedloop/cross_pairs" begin
    T, h, k_on, k_off, r, A = 0.01, 0.01 / 8, 50.0, 20.0, 0.01, 40.96
    k̃ = (1 - exp(-k_on * h)) / h
    for (case, DB) in ((1, 0.0), (2, 4.0))
        w = cross_world(30 + case, DB)
        num = 0.0; den = 0.0; nab = 0
        first_xy = Dict{Int,Tuple{Float64,Float64}}()
        bb = Set{Int}()          # B ids ever in a B-B pair (the non-keeper is placed d_dimer away)
        last_bound = Dict{Int,Tuple{Int,Float64,Float64}}()
        bad_links = 0; bad_still = 0
        for k in 1:2050
            SMLMSim.step!(w, (k - 1) * T, k * T)
            rows = frame_truth(w)
            id = byid(rows)
            if case == 1
                for r_ in rows
                    r_.pop == 2 || continue
                    r_.partner_pop == 2 && push!(bb, r_.id)
                    r_.partner_pop == 2 && push!(bb, r_.partner)
                    haskey(first_xy, r_.id) || (first_xy[r_.id] = (r_.x, r_.y))
                end
            end
            k > 50 || continue
            na = count(r_ -> r_.pop == 1 && r_.partner == 0, rows)
            nb = count(r_ -> r_.pop == 2 && r_.partner == 0, rows)
            nab_k = count(r_ -> r_.pop == 1 && r_.partner_pop == 2, rows)
            nab += nab_k
            num += k_off * nab_k
            den += k̃ * (π * r^2 / A) * na * nb
            for r_ in rows
                if r_.pop == 1 && r_.partner_pop == 2
                    o = id[r_.partner]
                    (o.partner == r_.id && o.partner_pop == 1) || (bad_links += 1)
                    if case == 1 && r_.bound == 1
                        prev = get(last_bound, r_.id, (0, 0.0, 0.0))
                        prev[1] == r_.partner && (r_.x, r_.y) != (prev[2], prev[3]) && (bad_still += 1)
                        last_bound[r_.id] = (r_.partner, r_.x, r_.y)
                        continue
                    end
                end
                delete!(last_bound, r_.id)
            end
        end
        nf = k_off * T * nab
        @test nf > 2000
        @test abs(num / den - 1) <= 3 * sqrt(2 / nf)
        @test bad_links == 0
        if case == 1
            w2 = cross_world(31, 0.0)   # re-run the first frame: B anchors keep their first-frame positions
            SMLMSim.step!(w2, 0.0, T)
            first = byid(frame_truth(w2))
            moved = 0
            ws = cross_world(31, 0.0)
            for k in 1:50
                SMLMSim.step!(ws, (k - 1) * T, k * T)
                for r_ in frame_truth(ws)
                    r_.pop == 2 && !(r_.id in bb) && haskey(first, r_.id) &&
                        (r_.x, r_.y) != (first[r_.id].x, first[r_.id].y) && (moved += 1)
                end
            end
            @test moved == 0
            @test bad_still == 0
        end
    end
end

@testset "closedloop/nonbinding" begin
    counts = Int[]
    for (binds_B, rep) in ((true, 1), (true, 2), (false, 1))
        w = cross_world(41, 0.0; k_on=Inf, binds_B=binds_B)
        nB = 0; bad_B = 0; nAA = 0
        for k in 1:500
            SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01)
            for r in frame_truth(w)
                if r.pop == 1
                    r.partner_pop == 2 && (nB += 1)
                    r.partner_pop == 1 && (nAA += 1)
                else
                    (r.partner != 0 && !binds_B) && (bad_B += 1)
                    (r.bound != 0 && !binds_B) && (bad_B += 1)
                end
            end
        end
        if binds_B
            @test nB > 1000
            push!(counts, nB)
        else
            @test nB == 0
            @test bad_B == 0
            @test nAA > 0
        end
    end
    @test counts[1] == counts[2]      # one seed, two runs
    dk = DimerKinetics(k_on=50.0, r_react=0.01, k_off=20.0, D_rot=1.0, d_dimer=0.05)
    mk(; kw...) = Population(density=1.0, fluor=one_state(1000.0), psf=GaussianPSF(0.1); kw...)
    build(p) = SimWorld(StableRNG(1), cam32(), [p]; n_sub=1, dimers=dk)
    @test_throws ArgumentError build(mk(multiplicity=10))
    @test build(mk(multiplicity=10, binds=false)) isa SimWorld
    @test build(mk(multiplicity=0)) isa SimWorld
    @test build(mk(multiplicity=1)) isa SimWorld
    @test SimWorld(StableRNG(1), cam32(), [mk(multiplicity=10)]; n_sub=1) isa SimWorld
    @test_throws ArgumentError SimWorld(StableRNG(1), IdealCamera(1:2, 1:2, 0.1), [mk()]; n_sub=1, margin=0.0,
                                        dimers=DimerKinetics(k_on=1.0, r_react=0.15, k_off=1.0, D_rot=0.0, d_dimer=0.0))
    @test_throws ArgumentError DimerKinetics(k_on=0.0, r_react=0.01, k_off=1.0, D_rot=0.0, d_dimer=0.0)
    @test_throws ArgumentError DimerKinetics(k_on=1.0, r_react=0.01, k_off=1.0, D_rot=0.0, d_dimer=0.0, D_dimer=:max)
    @test DimerKinetics(k_on=1, r_react=0.01, k_off=0, D_rot=0, d_dimer=0, D_dimer=1).D_dimer === 1.0
    cfg = DiffusionSMLMConfig(density=1.0, diff_dimer=0.07, k_off=0.3, r_react=0.02, d_dimer=0.04, diff_dimer_rot=0.6)
    d2 = DimerKinetics(cfg; k_on=9.0)
    @test (d2.k_on, d2.r_react, d2.k_off, d2.D_dimer, d2.D_rot, d2.d_dimer) == (9.0, 0.02, 0.3, 0.07, 0.6, 0.04)
end

@testset "closedloop/pair_departure" begin
    T = 0.01
    dk = DimerKinetics(k_on=Inf, r_react=0.05, k_off=5.0, D_rot=1.0, d_dimer=0.02)
    mover = Population(name=:mover, density=30.0, lifetime=0.05, mobility=[(1.0, 1.0)], fluor=one_state(1000.0),
                       psf=GaussianPSF(0.1))
    w = SimWorld(StableRNG(51), cam32(), [mover]; n_sub=4, margin=0.0, dimers=dk)
    counts = Float64[]
    bad = 0; ndep = 0
    for k in 1:3000
        SMLMSim.step!(w, (k - 1) * T, k * T)
        push!(counts, w.pops[1].n)
        rows = frame_truth(w)
        id = byid(rows)
        for r in rows
            (isnan(r.t_depart) || r.partner == 0) && continue
            ndep += 1
            id[r.partner].t_depart == r.t_depart || (bad += 1)
        end
    end
    @test ndep > 100
    @test bad == 0
    A = 3.2^2
    batch = [mean(counts[(b-1)*150+1:b*150]) for b in 1:20]
    @test abs(mean(batch) - 30.0 * A) <= 3 * std(batch) / sqrt(20)
    anchor = Population(name=:anchor, density=5.0, fluor=one_state(1000.0), psf=GaussianPSF(0.1))
    w2 = SimWorld(StableRNG(52), cam32(), [mover, anchor]; n_sub=4, margin=0.0, dimers=dk)
    n_anchor = Int[]
    bad2 = 0; nbound = 0
    for k in 1:1000
        SMLMSim.step!(w2, (k - 1) * T, k * T)
        rows = frame_truth(w2)
        push!(n_anchor, count(r -> r.pop == 2, rows))
        for r in rows
            r.pop == 1 && r.partner_pop == 2 && (nbound += 1)
            r.pop == 1 && r.partner_pop == 2 && !isnan(r.t_depart) && (bad2 += 1)
        end
    end
    @test allequal(n_anchor)
    @test nbound > 100
    @test bad2 == 0
end

@testset "closedloop/dimer_times" begin
    T = 0.01
    # (1) many forced pairs: formed at t0, break times exponential, links symmetric, bound exact
    k_off = 5.0
    dk = DimerKinetics(k_on=Inf, r_react=0.01, k_off=k_off, D_rot=0.0, d_dimer=0.005)
    pop = Population(density=0.0, fluor=one_state(1000.0), psf=GaussianPSF(0.05))
    w = SimWorld(StableRNG(61), cam32(), [pop]; n_sub=2, margin=0.0, dimers=dk)
    ps = w.pops[1]
    npair = 2000
    xs = Float64[]; ys = Float64[]
    for a in 1:64, b in 1:32
        length(xs) >= 2npair && break
        x = 0.05 * a; y = 0.1 * b
        append!(xs, (x, x + 0.004)); append!(ys, (y, y))
    end
    for _ in 1:2npair
        _add_emitter!(w, ps, 0.0)
    end
    place!(ps, xs, ys)
    tbreak = Dict{Int,Float64}()
    bad_sym = 0; bad_bound = 0; bad_form = 0
    for k in 1:300
        SMLMSim.step!(w, (k - 1) * T, k * T)
        rows = frame_truth(w)
        id = byid(rows)
        for r in rows
            k == 1 && r.t_form != 0.0 && (bad_form += 1)
            isnan(r.t_break) || (tbreak[r.id] = r.t_break)
            if r.partner != 0
                o = id[r.partner]
                (o.partner == r.id) || (bad_sym += 1)
            end
            tbk = get(tbreak, r.id, Inf)
            want = clamp(min(k * T, tbk) - (k - 1) * T, 0.0, T) / T
            abs(r.bound - want) <= 1e-12 || (bad_bound += 1)
        end
    end
    @test bad_sym == 0
    @test bad_bound == 0
    @test bad_form == 0
    tb = collect(values(tbreak))
    @test length(tb) == 2npair
    @test ks_exp(tb[1:2:end], 1 / k_off) * sqrt(length(tb) / 2) < 1.63
    @test tbreak[ps.id[1]] == tbreak[ps.id[2]]
    # (2) bleaching: the pair stays, dark, until its break; exact t_bleach of the first member
    dk2 = DimerKinetics(k_on=Inf, r_react=0.01, k_off=0.2, D_rot=0.0, d_dimer=0.005)
    pop2 = Population(density=0.0, fluor=one_state(1e4), psf=GaussianPSF(0.05))
    w2 = pair_world([pop2], dk2; seed=62)
    ps2 = w2.pops[1]
    ps2.budget[1] = 30.0        # bleaches at 3 ms, in frame 1
    ps2.budget[2] = 250.0       # bleaches at 25 ms, in frame 3
    id1, id2 = ps2.id[1], ps2.id[2]
    kb = 100
    bad = 0
    for k in 1:100
        SMLMSim.step!(w2, (k - 1) * T, k * T)
        rows = frame_truth(w2)
        length(rows) == 2 || (bad += 1)
        id = byid(rows)
        r1, r2 = id[id1], id[id2]
        k == 1 && abs(r1.t_bleach - 30.0 / 1e4) > 1e-12 && (bad += 1)
        k == 3 && abs(r2.t_bleach - 250.0 / 1e4) > 1e-12 && (bad += 1)
        k >= 1 && r1.m != 0 && (bad += 1)
        k >= 3 && r2.m != 0 && (bad += 1)
        any(r -> !isnan(r.t_depart), rows) && (bad += 1)
        if any(r -> !isnan(r.t_break), rows)
            kb = k
            break
        end
        (r1.partner == id2 && r2.partner == id1) || (bad += 1)
        want = 1.0     # formed at 0, no break yet
        (r1.bound == want && r2.bound == want) || (bad += 1)
    end
    @test bad == 0
    @test kb > 3
end

@testset "closedloop/unbind_no_reform" begin
    dk = DimerKinetics(k_on=Inf, r_react=0.01, k_off=20.0, D_rot=0.0, d_dimer=0.005)
    pop = Population(density=0.0, fluor=one_state(1000.0), psf=GaussianPSF(0.05))
    T = 0.01
    function run(w, x_first)
        nform = 0; nbreak = 0; ksplit = 0
        anchor_xy = nothing
        bad_late = 0; bad_sep = 0; bad_box = 0; bad_anchor = 0
        for k in 1:100
            SMLMSim.step!(w, (k - 1) * T, k * T)
            rows = sort(collect(frame_truth(w)); by=r -> r.id)
            any(r -> !isnan(r.t_form), rows) && (nform += 1)
            if any(r -> !isnan(r.t_break), rows)
                nbreak += 1
                ksplit = k
            end
            if k == 1
                @test all(r -> r.t_form == 0.0, rows)
                anchor_xy = (rows[1].x, rows[1].y)
            end
            abs(rows[1].x - anchor_xy[1]) > 1e-12 && (bad_anchor += 1)
            abs(rows[1].y - anchor_xy[2]) > 1e-12 && (bad_anchor += 1)
            if ksplit != 0 && k > ksplit
                all(r -> r.partner == 0 && r.bound == 0, rows) || (bad_late += 1)
                hypot(rows[1].x - rows[2].x, rows[1].y - rows[2].y) > 0.01 || (bad_sep += 1)
                all(r -> 0 <= r.x <= 3.2 && 0 <= r.y <= 3.2, rows) || (bad_box += 1)
            end
        end
        return nform, nbreak, ksplit, bad_late, bad_sep, bad_box, bad_anchor
    end
    nform, nbreak, ksplit, bad_late, bad_sep, bad_box, bad_anchor = run(pair_world([pop], dk; seed=71), nothing)
    @test nform == 1
    @test nbreak == 1
    @test 2 <= ksplit <= 50
    @test (bad_late, bad_sep, bad_box, bad_anchor) == (0, 0, 0, 0)
    nform, nbreak, ksplit, bad_late, bad_sep, bad_box, bad_anchor =
        run(pair_world([pop], dk; seed=71, xs=[0.002, 0.001], ys=[1.6, 1.6]), nothing)
    @test nform == 1
    @test nbreak == 1
    @test (bad_late, bad_sep, bad_box, bad_anchor) == (0, 0, 0, 0)
end

# analytic re-pair probability: the Gaussian mass of the disk |q| < r around the partner, seen from r0
function repair_p(r, r0, σ2)
    n = 400
    s = 0.0
    for a in 1:n, b in 1:n
        ρ = (a - 0.5) / n * r
        φ = (b - 0.5) / n * 2π
        qx = ρ * cos(φ)
        qy = ρ * sin(φ)
        s += exp(-((qx - r0)^2 + qy^2) / (2σ2)) / (2π * σ2) * ρ
    end
    return s * (r / n) * (2π / n)
end

@testset "closedloop/repair_rate" begin
    T = 0.01
    for (r, dd, ρ, nfr, seed) in ((0.03, 0.01, 0.5, 5000, 81), (0.3, 0.05, 0.02, 30000, 82))
        pa = Population(name=:A, density=ρ, mobility=[(1.0, 0.37)], fluor=one_state(1000.0), psf=GaussianPSF(0.1))
        pb = Population(name=:B, density=ρ, mobility=[(1.0, 0.0)], fluor=one_state(1000.0), psf=GaussianPSF(0.1))
        w = SimWorld(StableRNG(seed), IdealCamera(1:40, 1:40, 1.0), [pa, pb]; n_sub=1, margin=0.0,
                     boundary=:periodic,
                     dimers=DimerKinetics(k_on=Inf, r_react=r, k_off=10.0, D_rot=0.0, d_dimer=dd, D_dimer=:min))
        for k in 1:100
            SMLMSim.step!(w, (k - 1) * T, k * T)
        end
        nsplit = 0; nrep = 0
        pending = Tuple{Int,Int}[]
        for k in 101:100+nfr
            SMLMSim.step!(w, (k - 1) * T, k * T)
            rows = frame_truth(w)
            if !isempty(pending)
                id = byid(rows)
                tf = (k - 101 + 100) * T      # the split frame's end = this frame's start
                for (a, b) in pending
                    ra, rb = id[a], id[b]
                    ok = ra.t_form == tf && rb.t_form == tf &&
                         (ra.partner == b || (!isnan(ra.t_break) && ra.t_break == rb.t_break))
                    ok && (nrep += 1)
                end
                empty!(pending)
            end
            bt = Dict{Float64,Int}()
            for r_ in rows
                r_.pop == 2 && !isnan(r_.t_break) && (bt[r_.t_break] = r_.id)
            end
            for r_ in rows
                r_.pop == 1 && !isnan(r_.t_break) && haskey(bt, r_.t_break) &&
                    push!(pending, (r_.id, bt[r_.t_break]))
            end
            nsplit += length(pending)
        end
        p = repair_p(r, r * (1 + 1e-9), 2 * 0.37 * T)
        @test nsplit >= 2000
        @test abs(nrep / nsplit - p) <= 3 * sqrt(p * (1 - p) / nsplit)
    end
end

@testset "closedloop/partner_links" begin
    T = 0.01
    dk = DimerKinetics(k_on=Inf, r_react=0.05, k_off=5.0, D_rot=1.0, d_dimer=0.02)
    pa = Population(name=:A, density=20.0, lifetime=0.05, mobility=[(1.0, 1.0)], fluor=one_state(1000.0),
                    psf=GaussianPSF(0.1))
    pb = Population(name=:B, density=20.0, lifetime=0.2, budget=300.0, mobility=[(1.0, 0.5)],
                    fluor=one_state(1e5), psf=GaussianPSF(0.1))
    w = SimWorld(StableRNG(91), cam32(), [pa, pb]; n_sub=4, margin=0.0, dimers=dk)
    bad_ps = 0; bad_rows = 0; nlinked = 0
    for k in 1:300
        SMLMSim.step!(w, (k - 1) * T, k * T)
        for (kk, ps) in enumerate(w.pops), i in 1:ps.n
            ps.partner[i] == 0 && continue
            nlinked += 1
            kp, ip = Int(ps.partner_pop[i]), Int(ps.partner[i])
            ok = 1 <= kp <= 2 && 1 <= ip <= w.pops[kp].n &&
                 w.pops[kp].partner[ip] == i && w.pops[kp].partner_pop[ip] == kk &&
                 ps.partner_id[i] == w.pops[kp].id[ip]
            ok || (bad_ps += 1)
        end
        rows = frame_truth(w)
        id = byid(rows)
        for r in rows
            r.partner == 0 && continue
            (haskey(id, r.partner) && id[r.partner].partner == r.id) || (bad_rows += 1)
        end
    end
    @test nlinked > 100
    @test bad_ps == 0
    @test bad_rows == 0
end

function det_world(seed)
    flu = two_state(3000.0, 200.0, 20.0)
    tbl = StampTable(GaussianPSF(0.39), 0.1, range(-0.3, 1.0, length=8); radius=12)
    sig = Population(name=:sig, density=1.0, lifetime=0.1, fluor=flu, budget=100.0,
                     mobility=[(0.5, 0.3), (0.5, 0.05)], psf=GaussianPSF(0.13))
    oof = Population(name=:oof, layer=:oof, density=1.0, lifetime=0.05, fluor=flu, budget=200.0,
                     z=(0.5, 1.0), psf=tbl, binds=false)
    oof2 = Population(name=:oof2, layer=:oof, density=1.0, lifetime=0.05, fluor=flu, budget=200.0,
                      z=(-0.2, 0.2), psf=GaussianPSF(0.39), binds=false)
    anchor = Population(name=:anchor, density=1.0, fluor=flu, psf=GaussianPSF(0.13))
    bg = BackgroundModel(level=5000.0, jitter=0.05, contrast=0.3, feature_size=0.3,
                         correlation_time=0.05, illumination_width=2.0, stretch=0.1)
    dk = DimerKinetics(k_on=50.0, r_react=0.05, k_off=5.0, D_rot=1.0, d_dimer=0.02)
    return SimWorld(StableRNG(seed), cam32(), [sig, oof, oof2, anchor]; n_sub=4, background=bg, dimers=dk)
end

const DIMER_GOLDEN = 50337.671348101125

@testset "closedloop/dimer_determinism" begin
    a, b = det_world(21), det_world(21)
    T = 0.01
    same = true
    tot_bound = 0.0
    golden = 0.0
    for k in 1:50
        ea = copy(SMLMSim.step!(a, (k - 1) * T, k * T))
        eb = SMLMSim.step!(b, (k - 1) * T, k * T)
        same &= ea == eb && collect(frame_truth(a)) == collect(frame_truth(b))
        tot_bound += sum(r.bound for r in frame_truth(a))
        if k == 20
            golden = sum(ea) + sum(r.bound for r in frame_truth(a))
        end
    end
    @test same
    @test tot_bound > 0
    @test isapprox(golden, DIMER_GOLDEN; rtol=1e-12)
end

@testset "closedloop/zero_alloc_dimers" begin
    w = det_world(23)
    exc = PeriodicExc(0.0173)
    scam = SCMOSCamera(32, 32, 0.1, 1.6; offset=100.0, gain=2.0, qe=1.0)
    nrng = StableRNG(99)
    dst = zeros(32, 32)
    function count_allocs(k)
        return @allocated begin
            SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01, exc)
            SMLMSim.scmos_noise!(nrng, copyto!(dst, w.expected), scam)
        end
    end
    for k in 1:200
        SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01, exc)
    end
    count_allocs(201)
    total = 0
    for k in 202:301
        total += count_allocs(k)
    end
    @test total == 0
end

@testset "closedloop/dark_binders" begin
    T = 0.01
    dk = DimerKinetics(k_on=50.0, r_react=0.01, k_off=20.0, D_rot=1.0, d_dimer=0.005)
    pop = Population(density=10.0, mobility=[(1.0, 0.5)], fluor=one_state(1e4), budget=50.0, psf=GaussianPSF(0.05))
    function dark_run(dimers, nframes)
        w = SimWorld(StableRNG(101), cam32(), [pop]; n_sub=4, margin=0.0, dimers)
        frames = Vector{Vector{FrameTruth}}()
        for k in 1:nframes
            SMLMSim.step!(w, (k - 1) * T, k * T)
            push!(frames, collect(frame_truth(w)))
        end
        return frames
    end
    fr = dark_run(dk, 100)
    @test allequal(length.(fr))
    @test all(f -> all(r -> r.m == 0 && r.photons == 0, f), fr[10:end])
    a, b = byid(fr[50]), byid(fr[51])
    unb = [id for (id, r) in a if r.partner == 0 && r.bound == 0 && b[id].partner == 0]
    @test length(unb) > 50
    @test count(id -> (a[id].x, a[id].y) != (b[id].x, b[id].y), unb) > 0.9 * length(unb)
    @test any(r -> r.partner != 0, Iterators.flatten(fr[51:100]))
    fr0 = dark_run(nothing, 100)
    @test length(fr0[100]) < length(fr0[1])
    # (b) unlabeled molecule with a lit one
    lit = Population(name=:Y, density=0.0, fluor=one_state(1e4), psf=GaussianPSF(0.05))
    unl = Population(name=:X, density=0.0, multiplicity=0, fluor=one_state(1e4), psf=GaussianPSF(0.05))
    dkb = DimerKinetics(k_on=Inf, r_react=0.01, k_off=0.0, D_rot=0.0, d_dimer=0.005)
    w = pair_world([unl, lit], dkb; xs=[1.6, 1.604], ys=[1.6, 1.6])
    idx, idy = w.pops[1].id[1], w.pops[2].id[1]
    bad = 0
    for k in 1:20
        SMLMSim.step!(w, (k - 1) * T, k * T)
        id = byid(frame_truth(w))
        rx, ry = id[idx], id[idy]
        k == 1 && !(rx.t_form == 0.0 && ry.t_form == 0.0 && rx.partner == idy && ry.partner == idx) && (bad += 1)
        (rx.photons == 0 && rx.m == 0 && isnan(rx.t_bleach)) || (bad += 1)
        (ry.bound == 1 && isapprox(ry.photons, 1e4 * T; rtol=1e-12)) || (bad += 1)
    end
    @test bad == 0
    # (c) dark with dark
    w = pair_world([unl], dkb)
    bad = 0
    for k in 1:20
        SMLMSim.step!(w, (k - 1) * T, k * T)
        rows = frame_truth(w)
        (length(rows) == 2 && all(r -> r.bound == 1, rows)) || (bad += 1)
    end
    @test bad == 0
end

@testset "closedloop/visible_formation" begin
    T, γ = 0.01, 1e4
    dk = DimerKinetics(k_on=Inf, r_react=0.05, k_off=0.0, D_rot=0.0, d_dimer=0.02)
    fl = two_state(γ, 1e-6, 1e-6)
    function vis_world(xs, ys; pops=nothing, merge_radius=0.0, mult=(1, 1), dkin=dk)
        pops === nothing && (pops = [Population(name=:A, density=0.0, fluor=fl, psf=GaussianPSF(0.05)),
                                     Population(name=:B, density=0.0, fluor=fl, multiplicity=mult[2], psf=GaussianPSF(0.05))])
        w = SimWorld(StableRNG(7), cam32(), pops; n_sub=4, margin=0.0, dimers=dkin, merge_radius)
        for (k, ps) in enumerate(w.pops)
            _add_emitter!(w, ps, 0.0)
            place!(ps, xs[k:k], ys[k:k])
            ps.state[1] = 1
        end
        return w
    end
    frame!(w, k) = (SMLMSim.step!(w, (k - 1) * T, k * T); byid(frame_truth(w)))
    # (i) both emitting
    w = vis_world([1.6, 1.63], [1.6, 1.6])
    ia, ib = w.pops[1].id[1], w.pops[2].id[1]
    f = frame!(w, 1)
    for r in (f[ia], f[ib])
        @test r.t_form == 0.0
        @test r.vis_form
        @test isapprox(r.lit_bound, 1; atol=1e-12)
        @test r.vis_bound
    end
    # (ii) B in state 2
    w = vis_world([1.6, 1.63], [1.6, 1.6])
    w.pops[2].state[1] = 2
    ia, ib = w.pops[1].id[1], w.pops[2].id[1]
    f = frame!(w, 1)
    @test !f[ia].vis_form && !f[ib].vis_form
    @test f[ib].lit_bound == 0
    @test !f[ia].vis_bound && !f[ib].vis_bound
    @test isapprox(f[ia].lit_bound, 1; atol=1e-12)
    # (iii) B bleaches at 0.015 s
    w = vis_world([1.6, 1.63], [1.6, 1.6])
    w.pops[2].budget[1] = γ * 0.015
    ia, ib = w.pops[1].id[1], w.pops[2].id[1]
    f1 = frame!(w, 1)
    f2 = frame!(w, 2)
    @test isapprox(f2[ib].t_bleach, 0.015; atol=1e-12)
    @test f2[ib].m == 0
    # frame 2: B is lit while bound for 5 ms of the 10 ms exposure, A for all of it; both rows are a visible pair
    @test isapprox(f2[ib].lit_bound, 0.5; atol=1e-12) && isapprox(f2[ia].lit_bound, 1; atol=1e-12)
    @test f2[ia].vis_bound && f2[ib].vis_bound
    f3 = frame!(w, 3)
    @test f3[ib].m == 0 && f3[ib].lit_bound == 0
    @test f3[ia].partner == ib && f3[ib].partner == ia
    @test !f3[ia].vis_bound && !f3[ib].vis_bound
    # (iv) dark formation: B bleaches before contact
    w = vis_world([1.6, 2.1], [1.6, 1.6])
    w.pops[2].budget[1] = γ * 0.005
    ia, ib = w.pops[1].id[1], w.pops[2].id[1]
    f = frame!(w, 1)
    @test isapprox(f[ib].t_bleach, 0.005; atol=1e-12) && f[ib].m == 0
    @test all(isnan(r.t_form) for r in values(f))
    w.pops[2].x[1] = 1.63; w.pops[2].y[1] = 1.6
    f = frame!(w, 2)
    for r in (f[ia], f[ib])
        @test isapprox(r.t_form, 0.01; atol=1e-12)
        @test !r.vis_form
        @test !r.vis_bound
    end
    @test f[ia].partner == ib && f[ib].partner == ia
    @test f[ib].lit_bound == 0 && f[ib].m == 0
    @test isapprox(f[ia].lit_bound, 1; atol=1e-12)
    f = frame!(w, 3)
    for r in (f[ia], f[ib])
        @test isnan(r.t_form) && !r.vis_form && !r.vis_bound
    end
    @test f[ia].partner == ib && f[ib].partner == ia
    # (v) overlap
    pops3 = [Population(name=:A, density=0.0, fluor=fl, psf=GaussianPSF(0.05)), Population(name=:B, density=0.0, fluor=fl, psf=GaussianPSF(0.05)),
             Population(name=:C, density=0.0, fluor=fl, binds=false, psf=GaussianPSF(0.05))]
    w = vis_world([1.6, 1.63], [1.6, 1.6]; pops=pops3[1:2], merge_radius=0.25)
    f = frame!(w, 1)
    @test all(!r.overlap for r in values(f))
    w = vis_world([1.6, 1.63, 1.7], [1.6, 1.6, 1.6]; pops=pops3, merge_radius=0.25)
    f = frame!(w, 1)
    @test length(f) == 3 && all(r.overlap for r in values(f))
    # (vi) unlabeled binder
    w = vis_world([1.6, 1.63], [1.6, 1.6]; mult=(1, 0))
    ia, ib = w.pops[1].id[1], w.pops[2].id[1]
    f = frame!(w, 1)
    @test f[ia].t_form == 0.0 && f[ib].t_form == 0.0
    @test !f[ia].vis_form && !f[ib].vis_form
    @test f[ib].lit_bound == 0 && isnan(f[ib].t_bleach)
    @test !f[ia].vis_bound && !f[ib].vis_bound
    @test isapprox(f[ia].lit_bound, 1; atol=1e-12)
    # (vii) mid-frame formation (finite k_on): bound and lit from t_form to the end of the exposure
    w = vis_world([1.6, 1.63], [1.6, 1.6]; dkin=DimerKinetics(k_on=300.0, r_react=0.05, k_off=0.0, D_rot=0.0, d_dimer=0.02))
    ia, ib = w.pops[1].id[1], w.pops[2].id[1]
    k = 1
    f = frame!(w, k)
    while isnan(f[ia].t_form) && k < 10
        k += 1
        f = frame!(w, k)
    end
    tf, ta = f[ia].t_form, (k - 1) * T
    @test ta < tf < ta + T && f[ib].t_form == tf
    for r in (f[ia], f[ib])
        @test r.vis_form && r.vis_bound
        @test isapprox(r.bound, (ta + T - tf) / T; atol=1e-12) && isapprox(r.lit_bound, (ta + T - tf) / T; atol=1e-12)
    end
    # (viii) breakup at 0.0137 s, inside frame 2: bound and lit for 3.7 ms of it, so a visible pair in it, while the
    # row's partner is the one at the end of the exposure (0); unbound in frame 3
    w = vis_world([1.6, 1.63], [1.6, 1.6])
    ia, ib = w.pops[1].id[1], w.pops[2].id[1]
    frame!(w, 1)
    w.pops[1].t_break_due[1] = 0.0137; w.pops[2].t_break_due[1] = 0.0137
    f = frame!(w, 2)
    for r in (f[ia], f[ib])
        @test isapprox(r.t_break, 0.0137; atol=1e-12)
        @test isapprox(r.bound, 0.37; atol=1e-9) && isapprox(r.lit_bound, 0.37; atol=1e-9)
        @test r.partner == 0 && r.vis_bound && !r.vis_form
    end
    f = frame!(w, 3)
    for r in (f[ia], f[ib])
        @test r.partner == 0 && r.lit_bound == 0 && !r.vis_bound
    end
    # (ix) B blinks off for frame 2 and on again for frame 3: not a visible pair in frame 2, again in frame 3;
    # vis_form only in the formation frame
    w = vis_world([1.6, 1.63], [1.6, 1.6])
    ia, ib = w.pops[1].id[1], w.pops[2].id[1]
    f = frame!(w, 1)
    @test f[ia].vis_bound && f[ib].vis_bound && f[ia].vis_form
    w.pops[2].state[1] = 2
    f = frame!(w, 2)
    @test f[ib].lit_bound == 0 && isapprox(f[ia].lit_bound, 1; atol=1e-12)
    @test !f[ia].vis_bound && !f[ib].vis_bound
    w.pops[2].state[1] = 1
    f = frame!(w, 3)
    @test isapprox(f[ib].lit_bound, 1; atol=1e-12) && f[ia].vis_bound && f[ib].vis_bound
    @test !f[ia].vis_form && !f[ib].vis_form
    # (x) the breakup frame of (viii) with merge_radius 0.25: a pair bound during the exposure is not an overlap
    # although it split in it; unbound through frame 3, the same two rows are
    w = vis_world([1.6, 1.63], [1.6, 1.6]; merge_radius=0.25)
    ia, ib = w.pops[1].id[1], w.pops[2].id[1]
    frame!(w, 1)
    w.pops[1].t_break_due[1] = 0.0137; w.pops[2].t_break_due[1] = 0.0137
    f = frame!(w, 2)
    for r in (f[ia], f[ib])
        @test r.partner == 0 && r.vis_bound && !r.overlap
    end
    f = frame!(w, 3)
    for r in (f[ia], f[ib])
        @test r.bound == 0 && r.overlap
    end
    # (xi) formed at 0 and broken inside frame 1 (k_off 1000/s): a visible pair in frame 1, partner 0 at its end
    w = vis_world([1.6, 1.63], [1.6, 1.6]; dkin=DimerKinetics(k_on=Inf, r_react=0.05, k_off=1000.0, D_rot=0.0, d_dimer=0.02))
    ia, ib = w.pops[1].id[1], w.pops[2].id[1]
    f = frame!(w, 1)
    for r in (f[ia], f[ib])
        @test r.t_form == 0.0 && 0 < r.t_break < T && r.partner == 0
        @test r.vis_form && r.vis_bound && isapprox(r.bound, r.t_break / T; atol=1e-12)
    end
    # (xii) a break due exactly at the start of frame 2: no bound time in frame 2, so not a visible pair there, and
    # with merge_radius 0.25 the two rows overlap
    w = vis_world([1.6, 1.63], [1.6, 1.6]; merge_radius=0.25)
    ia, ib = w.pops[1].id[1], w.pops[2].id[1]
    f = frame!(w, 1)
    @test f[ia].vis_bound && f[ib].vis_bound && !f[ia].overlap
    w.pops[1].t_break_due[1] = T; w.pops[2].t_break_due[1] = T
    f = frame!(w, 2)
    for r in (f[ia], f[ib])
        @test r.bound == 0 && r.lit_bound == 0 && !r.vis_bound && r.overlap
    end
end

# the design section 7a configuration itself (`bench_setup`: 256 x 256 sCMOS, about 250 in-focus binders with
# blinking, bleaching, births and dimers, about 20 OOF blobs, about 30 haze blobs, the Cell9-like background,
# merge_radius, a SpotExcitation with a switching spot): the measured steps 2.01-3.01 s include both switches
@testset "closedloop/zero_alloc_full" begin
    w, scam, exc = bench_setup(31)
    nrng = StableRNG(98)
    dst = zeros(BENCH_N, BENCH_N)
    function count_allocs(k)
        return @allocated begin
            SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01, exc)
            s = 0.0
            for r in frame_truth(w)
                s += r.photons + r.bound + (r.overlap ? 1.0 : 0.0)
            end
            SMLMSim.scmos_noise!(nrng, copyto!(dst, w.expected), scam)
            s
        end
    end
    for k in 1:200
        SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01, exc)
    end
    count_allocs(201)
    total = 0
    nmol = 0
    noof = 0
    npair = 0
    for k in 202:301
        total += count_allocs(k)
        for r in frame_truth(w)
            r.photons > 0 && w.pops[r.pop].p.name === :mol && (nmol += 1)
            r.photons > 0 && w.pops[r.pop].p.name === :oof && (noof += 1)
            r.partner != 0 && (npair += 1)
        end
    end
    @test total == 0
    # the workload: section 7a's render budget has every in-focus emitter and blob rendering
    @test nmol / 100 >= 250 && noof / 100 >= 20 && npair > 0
end

@testset "closedloop/params_dict" begin
    allowed(v) = v isa Union{String,Bool,Int,Float64,Vector{Float64}}
    haskey_prefix(d, k) = haskey(d, k) || any(s -> startswith(s, k * "."), keys(d))
    w = det_world(21)
    d = SMLMSim.params_dict(w)
    @test d isa Dict{String,Any}
    @test all(allowed, values(d))
    @test d["smlmsim.version"] == string(pkgversion(SMLMSim))
    @test d["world.n_sub"] == 4 && d["world.boundary"] == "reflecting" && d["world.merge_radius"] == 0.0
    @test d["world.box_um"] == collect(Float64, w.box)
    @test d["camera.type"] == "IdealCamera" && d["camera.nx"] == 32 && d["camera.ny"] == 32
    @test d["camera.pixel_size_um"] == 0.1
    for k in eachindex(w.pops), f in fieldnames(Population)
        @test haskey_prefix(d, "pop$k.$f")
    end
    @test d["pop1.binds"] === true && d["pop2.binds"] === false
    @test d["pop1.mobility.fraction"] == [0.5, 0.5] && d["pop1.mobility.D"] == [0.3, 0.05]
    # design section 6.6: the rate matrix as `.q` (row-major) with `.q.n`
    @test d["pop1.q"] == [-200.0, 200.0, 20.0, -20.0] && d["pop1.q.n"] == 2
    @test !haskey(d, "pop1.fluor.q") && !haskey(d, "pop1.fluor.q.n")
    @test d["pop1.fluor.gamma"] == 3000.0
    @test d["pop1.psf.sigma_um"] == 0.13
    @test d["pop2.psf.stamp.radius"] == 12 && d["pop2.psf.stamp.z_min"] ≈ -0.3 && d["pop2.psf.stamp.z_max"] ≈ 1.0
    @test d["pop2.psf.stamp.oversample"] == 4 && d["pop2.psf.stamp.z_step"] ≈ 1.3 / 7
    for k in ("k_on", "r_react", "k_off", "D_rot", "d_dimer")
        @test d["dimers.$k"] isa Float64
    end
    @test d["dimers.D_dimer"] == "min"
    for f in fieldnames(BackgroundModel)
        @test haskey(d, "bg.$f")
    end
    @test d["bg.level"] == 5000.0 && d["bg.illumination_width"] == 2.0
    # no dimers, no background, a fixed D_dimer, a distribution level and per-pixel sCMOS maps
    pop = Population(density=1.0, fluor=one_state(1e3), psf=GaussianPSF(0.1))
    scam = SCMOSCamera(32, 32, 0.1, fill(1.6, 32, 32); offset=fill(100.0, 32, 32), gain=2.0, qe=0.9)
    w0 = SimWorld(StableRNG(1), scam, [pop]; n_sub=2)
    d0 = SMLMSim.params_dict(w0)
    @test all(allowed, values(d0))
    @test !any(k -> startswith(k, "dimers.") || startswith(k, "bg."), keys(d0))
    @test d0["camera.type"] == "SCMOSCamera" && d0["camera.readnoise.mean"] ≈ 1.6 && d0["camera.offset.mean"] ≈ 100.0
    @test d0["camera.gain"] == 2.0 && d0["camera.qe"] == 0.9
    wd = SimWorld(StableRNG(1), cam32(), [pop]; n_sub=2, boundary=:periodic, merge_radius=0.25,
                  dimers=DimerKinetics(k_on=Inf, r_react=0.05, k_off=1.0, D_rot=0.0, d_dimer=0.02, D_dimer=0.5),
                  background=BackgroundModel(level=Uniform(1.0, 2.0)))
    dd = SMLMSim.params_dict(wd)
    @test all(allowed, values(dd))
    @test dd["dimers.D_dimer"] === 0.5 && dd["dimers.k_on"] == Inf && dd["world.boundary"] == "periodic"
    @test dd["world.merge_radius"] == 0.25
    @test dd["bg.level"] isa String && occursin("Uniform", dd["bg.level"])
end

# the 5 ms budget of design section 7a, measured only on request (kitt): SMLMSIM_BENCH=1
if get(ENV, "SMLMSIM_BENCH", "") == "1"
    @testset "closedloop/benchmark" begin
        r = bench_run()
        bench_report(r)
        @test r.both <= 5e-3
        @test r.alloc_step == 0 && r.alloc_sweep == 0 && r.alloc_noise == 0
        @test r.emitting_mol >= 250 && r.emitting_oof >= 20
    end
end

# ---- 0.8.0 fix (Codex r1 on lane A): frame truth ----

struct SwitchExc
    ts::Float64
    declared::Bool
end
(e::SwitchExc)(x::Float64, y::Float64, z::Float64, t::Float64) = t < e.ts ? 1.0 : 3.0
SMLMSim.next_switch(e::SwitchExc, t::Float64) = e.declared && t < e.ts ? e.ts : Inf

# a world of one immobile one_state(1000.0) emitter at (1.6, 1.6); `kw` goes to the Population
function ev_world(; n_sub=1, boundary=:reflecting, fluor=one_state(1000.0), t0=0.0, kw...)
    pop = Population(density=0.0, fluor=fluor, psf=GaussianPSF(0.05); kw...)
    w = SimWorld(StableRNG(1), cam32(), [pop]; n_sub, margin=0.0, boundary, t0)
    ps = w.pops[1]
    _add_emitter!(w, ps, t0)
    place!(ps, [1.6], [1.6])
    return w, ps
end

@testset "closedloop/boundary_events" begin
    # an event exactly at t_b belongs to the next exposure
    # (a) a departure at t_b
    w, ps = ev_world(n_sub=1)
    ps.t_depart[1] = 0.01
    SMLMSim.step!(w, 0.0, 0.01)
    r = only(frame_truth(w))
    @test isnan(r.t_depart)
    @test r.photons ≈ 10 rtol = 1e-12
    @test ps.n == 1
    SMLMSim.step!(w, 0.01, 0.02)
    r = only(frame_truth(w))
    @test r.t_depart == 0.01
    @test r.photons == 0
    @test ps.n == 0
    # (b) the same at the end of the 8th sub-step
    w, ps = ev_world(n_sub=8)
    ps.t_depart[1] = 0.02
    SMLMSim.step!(w, 0.0, 0.01)
    SMLMSim.step!(w, 0.01, 0.02)
    @test isnan(only(frame_truth(w)).t_depart)
    SMLMSim.step!(w, 0.02, 0.03)
    @test only(frame_truth(w)).t_depart == 0.02
    # (c) a bleach at t_b
    w, ps = ev_world(n_sub=1, budget=1e6)
    ps.budget[1] = 10.0
    SMLMSim.step!(w, 0.0, 0.01)
    r = only(frame_truth(w))
    @test r.photons ≈ 10 rtol = 1e-12
    @test isnan(r.t_bleach)
    @test r.m == 1
    SMLMSim.step!(w, 0.01, 0.02)
    r = only(frame_truth(w))
    @test r.t_bleach == 0.01
    @test r.m == 0
    @test r.photons == 0
    # (d) a birth at t_b
    w, ps = ev_world(n_sub=8, lifetime=1e3, birth_rate=1e-9)
    ps.n = 0
    SMLMSim.step!(w, 0.0, 0.01)
    ps.t_next_birth = 0.02
    SMLMSim.step!(w, 0.01, 0.02)
    @test isempty(frame_truth(w))
    SMLMSim.step!(w, 0.02, 0.03)
    @test only(frame_truth(w)).t_birth == 0.02
    # (e) guard: a departure at the end of a gap
    w, ps = ev_world(n_sub=4)
    ps.t_depart[1] = 0.05
    SMLMSim.step!(w, 0.0, 0.01)
    SMLMSim.step!(w, 0.05, 0.06)
    r = only(frame_truth(w))
    @test r.t_depart == 0.05
    @test r.photons == 0
end

@testset "closedloop/fractions_in_range" begin
    w, ps = ev_world(n_sub=10)
    SMLMSim.step!(w, 0.0, 0.01)
    r = only(frame_truth(w))
    @test 0 <= r.lit <= 1
    @test r.lit ≈ 1 rtol = 1e-12
    pop = Population(density=0.0, fluor=one_state(1000.0), psf=GaussianPSF(0.05))
    dk = DimerKinetics(k_on=Inf, r_react=0.01, k_off=0.0, D_rot=0.0, d_dimer=0.005)
    w = pair_world([pop], dk; n_sub=10)
    for k in 1:5
        SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01)
        for r in frame_truth(w)
            @test 0 <= r.bound <= 1
            @test 0 <= r.lit_bound <= 1
            k >= 2 && @test r.bound ≈ 1 rtol = 1e-12
        end
    end
end

@testset "closedloop/truth_integration" begin
    exc = UniformExcitation()
    # drive two recorded sub-steps by hand, moving the emitter between them
    function two_substeps(w, ps, x2, tm, exc)
        w.t_a, w.t_b = 0.0, 0.01
        _begin_truth!(w)
        _substep!(w, 0.0, tm, exc, true)
        ps.x[1] = x2
        _substep!(w, tm, 0.01, exc, true)
        _finish_truth!(w)
        return only(frame_truth(w))
    end
    # (a) photon-weighted centroid: weights 1:3
    w, ps = ev_world()
    ps.x[1] = 1.0
    r = two_substeps(w, ps, 2.0, 0.0025, exc)
    @test r.x ≈ 1.75 atol = 1e-12
    @test r.photons ≈ 10 atol = 1e-12
    # (b) no photons: the presence-weighted centroid
    w, ps = ev_world(fluor=GenericFluor(; γ=1000.0, q=[-1e-12 1e-12; 1e-12 -1e-12]))
    ps.state[1] = 2
    ps.x[1] = 1.0
    r = two_substeps(w, ps, 2.0, 0.005, exc)
    @test r.photons == 0
    @test r.x ≈ 1.5 atol = 1e-12
    # (c) the centroid across a periodic wall
    w, ps = ev_world(boundary=:periodic)
    ps.x[1] = 3.19
    r = two_substeps(w, ps, 0.01, 0.0025, exc)
    @test r.x ≈ 0.005 atol = 1e-12
    # (d) an emitter born and departed inside the exposure, across a declared switch
    w, ps = ev_world()
    ps.n = 0
    j = _add_emitter!(w, ps, 0.003)
    place!(ps, [1.6], [1.6])
    ps.t_depart[j] = 0.007
    w.t_a, w.t_b = 0.0, 0.01
    w.n_truth = 0
    e, alive, departed = _advance_acc!(w, ps, j, 0.0, 0.01, 0.003, SwitchExc(0.004, true))
    _write_row!(w, 1, ps, j, departed)
    r = only(frame_truth(w))
    @test departed
    @test r.photons ≈ 10 rtol = 1e-12
    @test r.lit ≈ 0.4 rtol = 1e-12
    @test r.excitation ≈ 2.5 rtol = 1e-12
    @test r.t_birth == 0.003
    @test r.t_depart == 0.007
    # (e) the accumulators reset after a gap
    w, ps = ev_world()
    SMLMSim.step!(w, 0.0, 0.01)
    SMLMSim.step!(w, 0.05, 0.06)
    r = only(frame_truth(w))
    @test r.photons ≈ 10 rtol = 1e-12
    @test r.lit ≈ 1 rtol = 1e-12
end

@testset "closedloop/bleach_in_frame" begin
    # a bleach inside the exposure: photons > 0 with m = 0 (the identity photons = m γ I lit T needs constant m)
    w, ps = ev_world(n_sub=4, budget=1e6)
    ps.budget[1] = 4.0
    SMLMSim.step!(w, 0.0, 0.01)
    r = only(frame_truth(w))
    @test r.photons ≈ 4 rtol = 1e-12
    @test r.m == 0
    @test r.lit ≈ 0.4 rtol = 1e-12
    @test r.t_bleach ≈ 0.004 rtol = 1e-12
end

# ---- 0.8.0 fix (Codex r1 on lane B): pairs ----

# one binding one_state(1000.0) population of immobile emitters, for the pair tests below
bpop(; kw...) = Population(density=0.0, fluor=one_state(1000.0), psf=GaussianPSF(0.05); kw...)

step_alloc(w, k) = @allocated SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01)

@testset "closedloop/pair_swap_order" begin
    # a swap-remove inside a bound pair's population must not exchange the members
    dk = DimerKinetics(k_on=Inf, r_react=0.02, k_off=0.0, D_rot=0.0, d_dimer=0.01, D_dimer=0.0)
    w = SimWorld(StableRNG(1), cam32(), [bpop()]; n_sub=2, margin=0.0, dimers=dk)
    ps = w.pops[1]
    for _ in 1:3
        _add_emitter!(w, ps, 0.0)
    end
    place!(ps, [2.5, 0.997, 1.007], [2.5, 1.6, 1.6])
    ps.t_depart[1] = 0.003
    SMLMSim.step!(w, 0.0, 0.01)
    SMLMSim.step!(w, 0.01, 0.02)
    x = Dict(ps.id[i] => ps.x[i] for i in 1:ps.n)
    @test x[2] ≈ 0.997 atol = 1e-12
    @test x[3] ≈ 1.007 atol = 1e-12
end

@testset "closedloop/dark_retention_bounds" begin
    dk = DimerKinetics(k_on=Inf, r_react=0.01, k_off=1.0, D_rot=0.0, d_dimer=0.005)
    g(; kw...) = Population(density=0.0, lifetime=Inf, budget=1.0, birth_rate=1000.0, fluor=one_state(1000.0),
                            psf=GaussianPSF(0.05); kw...)
    mkw(pops; kw...) = SimWorld(StableRNG(4), cam32(), pops; n_sub=1, margin=0.0, kw...)
    # (a) a binder with births and an infinite lifetime keeps its bleached molecules: no bound on the count
    @test_throws ArgumentError mkw([g()]; dimers=dk)
    @test mkw([g()]) isa SimWorld
    @test mkw([g(binds=false)]; dimers=dk) isa SimWorld
    # (b) births with an infinite lifetime need a molecule that can bleach
    @test_throws ArgumentError Population(density=0.0, birth_rate=10.0, budget=100.0, multiplicity=0,
                                          fluor=one_state(1000.0), psf=GaussianPSF(0.05))
    @test_throws ArgumentError Population(density=0.0, birth_rate=10.0, budget=100.0, multiplicity=1,
                                          fluor=one_state(0.0), psf=GaussianPSF(0.05))
    # (c) a non-binder is removed at bleaching even with dimers: the count stays bounded and steps allocate nothing
    w = mkw([g(binds=false)]; dimers=dk)
    total = 0
    for k in 1:300
        a = step_alloc(w, k)
        k > 200 && (total += a)
    end
    @test w.pops[1].n < 200
    @test total == 0
    # (d) a bleached binder stays, a bleached non-binder leaves
    pb = bpop(name=:b, budget=1e6)
    pn = bpop(name=:nb, budget=1e6, binds=false)
    w = mkw([pb, pn]; dimers=dk)
    for (k, ps) in enumerate(w.pops)
        _add_emitter!(w, ps, 0.0)
        place!(ps, [1.0 + 2.0 * (k - 1)], [1.6])
        ps.budget[1] = 4.0
    end
    SMLMSim.step!(w, 0.0, 0.01)
    rows = collect(frame_truth(w))
    @test length(rows) == 2
    @test all(r -> isapprox(r.t_bleach, 0.004; rtol=1e-9), rows)
    SMLMSim.step!(w, 0.01, 0.02)
    r = only(frame_truth(w))
    @test r.pop == 1 && r.m == 0 && r.photons == 0
end

@testset "closedloop/pair_reflection" begin
    # (a) a moving pair's centre reflects, a placed pair's is clamped
    w = SimWorld(StableRNG(1), cam32(), [bpop()]; n_sub=1, margin=0.0)
    ps = w.pops[1]
    for _ in 1:2
        _add_emitter!(w, ps, 0.0)
    end
    _put_pair!(w, ps, 1, ps, 2, 0.09, 1.6, 1.0, 0.0, 0.2, true)
    @test ps.x[1] ≈ 0.01 atol = 1e-12
    @test ps.x[2] ≈ 0.21 atol = 1e-12
    _put_pair!(w, ps, 1, ps, 2, 0.09, 1.6, 1.0, 0.0, 0.2, false)
    @test ps.x[1] ≈ 0.0 atol = 1e-12
    @test ps.x[2] ≈ 0.2 atol = 1e-12
    # (b) no bound member piles up on a wall
    pop = Population(density=20.0, mobility=[(1.0, 1.0)], fluor=one_state(1000.0), psf=GaussianPSF(0.05))
    dk = DimerKinetics(k_on=Inf, r_react=0.05, k_off=0.0, D_rot=1.0, d_dimer=0.04, D_dimer=1.0)
    w = SimWorld(StableRNG(2), cam32(), [pop]; n_sub=4, margin=0.0, dimers=dk)
    ps = w.pops[1]
    on_wall = 0
    outside = 0
    nbound = 0
    for k in 1:100
        SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01)
        for i in 1:ps.n
            outside += !(0.0 <= ps.x[i] <= 3.2 && 0.0 <= ps.y[i] <= 3.2)
            if ps.partner[i] != 0
                nbound += 1
                on_wall += ps.x[i] == 0.0 || ps.y[i] == 0.0
            end
        end
    end
    @test nbound > 1000
    @test on_wall == 0
    @test outside == 0
end

@testset "closedloop/dimer_validation" begin
    mk(; kw...) = DimerKinetics(; k_on=Inf, r_react=0.01, k_off=0.0, D_rot=0.0, d_dimer=0.005, kw...)
    @test mk() isa DimerKinetics
    for kw in ((D_rot=Inf,), (D_rot=NaN,), (D_dimer=Inf,), (D_dimer=NaN,), (d_dimer=Inf,))
        @test_throws ArgumentError mk(; kw...)
    end
    # a number that overflows or underflows Float64 is checked after the conversion
    for kw in ((D_rot=big"1e400",), (D_dimer=big"1e400",), (d_dimer=big"1e400",), (k_off=big"1e400",),
               (r_react=big"1e400",), (r_react=big"1e-400",))
        @test_throws ArgumentError mk(; kw...)
    end
    @test mk(; k_on=big"1e400").k_on == Inf
end

@testset "closedloop/split_margin" begin
    # coincident partners split to just over r_react, a distance that survives rounding: no re-formation
    for D_dimer in (:min, 0.0)
        dk = DimerKinetics(k_on=Inf, r_react=1e-8, k_off=1e4, D_rot=0.0, d_dimer=0.0, D_dimer=D_dimer)
        w = SimWorld(StableRNG(1), cam32(), [bpop()]; n_sub=1, margin=0.0, dimers=dk)
        ps = w.pops[1]
        for _ in 1:2
            _add_emitter!(w, ps, 0.0)
        end
        place!(ps, [1.6, 1.6], [1.6, 1.6])
        for k in 1:10
            SMLMSim.step!(w, (k - 1) * 0.01, k * 0.01)
            if k == 1
                @test _norm((_disp(ps.x[1], ps.x[2], nothing), _disp(ps.y[1], ps.y[2], nothing))) >= 1e-8
            else
                @test all(r -> r.partner == 0 && isnan(r.t_form), frame_truth(w))
            end
        end
    end
end

@testset "closedloop/periodic_axis" begin
    # the periodic displacement keeps the subtraction's rounding error: 1e-17 and 3.2 are 1e-17 apart
    dk = DimerKinetics(k_on=Inf, r_react=0.01, k_off=0.0, D_rot=0.0, d_dimer=0.005)
    w = pair_world([bpop()], dk; boundary=:periodic, xs=[1e-17, 3.2], ys=[1.6, 1.6])
    SMLMSim.step!(w, 0.0, 0.01)
    ps = w.pops[1]
    @test ps.x[1] == 1e-17
    @test ps.x[2] ≈ 3.195 atol = 1e-12
end

@testset "closedloop/tiny_scale_contact" begin
    cam = IdealCamera(1:32, 1:32, 1e-200)
    pop = Population(density=0.0, fluor=one_state(1000.0), psf=GaussianPSF(1e-201))
    dk = DimerKinetics(k_on=Inf, r_react=1e-200, k_off=0.0, D_rot=0.0, d_dimer=0.0)
    w = SimWorld(StableRNG(1), cam, [pop]; n_sub=1, margin=0.0, dimers=dk)
    ps = w.pops[1]
    for _ in 1:2
        _add_emitter!(w, ps, 0.0)
    end
    place!(ps, [1e-199, 1.05e-199], [1.6e-199, 1.6e-199])
    SMLMSim.step!(w, 0.0, 0.01)
    rows = frame_truth(w)
    @test length(rows) == 2
    @test rows[1].partner == rows[2].id && rows[2].partner == rows[1].id
    @test all(r -> r.t_form == 0.0, rows)
end

@testset "closedloop/cell_count" begin
    pop = Population(density=1.0, fluor=one_state(1000.0), psf=GaussianPSF(0.1))
    dk = DimerKinetics(k_on=Inf, r_react=1e-20, k_off=0.0, D_rot=0.0, d_dimer=0.0)
    w = SimWorld(StableRNG(1), cam32(), [pop]; n_sub=1, margin=0.0, dimers=dk)
    @test w.ncx == w.ncy == 128
    @test SMLMSim.step!(w, 0.0, 0.01) isa Matrix{Float64}
end

@testset "closedloop/formation_timing" begin
    # a finite k_on forms at the exact time E/k_on (E the first draw of the step) when it falls inside the sub-step
    mk(k_on) = pair_world([bpop()], DimerKinetics(k_on=k_on, r_react=0.01, k_off=0.0, D_rot=0.0, d_dimer=0.005); n_sub=1)
    ζ = randexp(copy(mk(50.0).rng))
    w = mk(ζ / 0.006)
    SMLMSim.step!(w, 0.0, 0.01)
    rows = frame_truth(w)
    @test length(rows) == 2
    @test all(r -> isapprox(r.t_form, 0.006; rtol=1e-12) && isapprox(r.bound, 0.4; rtol=1e-9), rows)
    @test rows[1].partner == rows[2].id && rows[2].partner == rows[1].id
    # a departure before the formation time prevents the pair
    w2 = mk(ζ / 0.006)
    w2.pops[1].t_depart[2] = 0.005
    SMLMSim.step!(w2, 0.0, 0.01)
    rows = frame_truth(w2)
    @test all(r -> isnan(r.t_form), rows)
    @test count(r -> r.t_depart == 0.005, rows) == 1
end

# ---- 0.8.0 fix (Codex r2 on lane A): an event belongs to the sub-step that contains its absolute time ----

# relative intensity 1 before `ts`, `hi` from `ts` on, a declared switch
struct JumpExc
    ts::Float64
    hi::Float64
end
(e::JumpExc)(x::Float64, y::Float64, z::Float64, t::Float64) = t < e.ts ? 1.0 : e.hi
SMLMSim.next_switch(e::JumpExc, t::Float64) = t < e.ts ? e.ts : Inf

# The pairs (b, time) of a monotone time function `f` on the floats around b0: when some b gives exactly `te`, the
# smallest such b; otherwise the largest b below and the smallest above, each with the time it gives, the nearest
# attainable on each side of an unattainable target.
function attain(f, b0, te; increasing=true)
    up(b) = increasing ? nextfloat(b) : prevfloat(b)
    down(b) = increasing ? prevfloat(b) : nextfloat(b)
    b = b0
    while f(b) < te
        b = up(b)
    end
    while f(b) >= te
        b = down(b)
    end
    hi = up(b)
    f(hi) == te && return [(hi, te)]
    return [(b, f(b)), (hi, f(hi))]
end

# Every pair (b, f(b)) with f(b) == te among the 8 floats either side of the smallest such b (the parameter values
# that map to one target time leave residuals of both signs); the nearest pairs around te when none attains it.
function attain_all(f, b0, te; increasing=true, ulps=8)
    base = attain(f, b0, te; increasing)
    cand = [base[1][1]]
    lo = hi = cand[1]
    for _ in 1:ulps
        lo = prevfloat(lo)
        hi = nextfloat(hi)
        pushfirst!(cand, lo)
        push!(cand, hi)
    end
    exact = [(b, f(b)) for b in cand if f(b) == te]
    return isempty(exact) ? base : exact
end

# the rows of the exposure so far, including the emitters still present, leaving the truth buffer as it was
function truth_snapshot(w)
    n = w.n_truth
    _finish_truth!(w)
    rows = copy(frame_truth(w))
    w.n_truth = n
    return rows
end

# the times at which the rows record an event of `kind` (a pair's two rows record one event)
boundary_events(kind, rows) =
    kind === :departure ? [r.t_depart for r in rows if !isnan(r.t_depart)] :
    kind === :birth ? [r.t_birth for r in rows if !isnan(r.t_birth)] :
    kind === :bleach || kind === :pending_bleach ? [r.t_bleach for r in rows if !isnan(r.t_bleach)] :
    kind === :break ? unique([r.t_break for r in rows if !isnan(r.t_break)]) :
    kind === :formation ? unique([r.t_form for r in rows if !isnan(r.t_form)]) : Float64[]

# Two consecutive windows A = [a0, a1) and B = [a1, b1) of one emitter world, as sub-steps driven by hand inside the
# frame [0, 0.01) (:interior) or as the frames A and B of step! (:frame), after one warm-up window ending at a0
# (none with warm = false). `setup!(w, ps)` places the event after the warm-up. Returns the events recorded in A
# and over A and B, the state of the emitter after each window and the photons of each.
function boundary_run(mkworld, setup!, exc, kind, mode, win; warm=true, run_b=true)
    a0, a1, b1 = win
    w, ps = mkworld()
    if mode === :interior
        w.t_a, w.t_b = 0.0, 0.01
        _begin_truth!(w)
        warm && _substep!(w, 0.0, a0, exc, true)
    else
        warm && SMLMSim.step!(w, a0 - 0.01, a0, exc)
    end
    setup!(w, ps)
    advance(t0, t1) = if mode === :interior
        _substep!(w, t0, t1, exc, true)
        truth_snapshot(w)
    else
        SMLMSim.step!(w, t0, t1, exc)
        copy(frame_truth(w))
    end
    ph0 = mode === :interior && ps.n > 0 ? ps.sph[1] : 0.0
    rowsA = advance(a0, a1)
    evA = boundary_events(kind, rowsA)
    stA = ps.n > 0 ? Int(ps.state[1]) : 0
    phA = ps.n > 0 ? ps.sph[1] - ph0 : NaN
    mA = isempty(rowsA) ? -1 : Int(rowsA[1].m)
    run_b || return (; evA, stA, phA, mA)
    ph1 = mode === :interior && ps.n > 0 ? ps.sph[1] : 0.0
    rowsB = advance(a1, b1)
    evT = mode === :interior ? boundary_events(kind, rowsB) : vcat(evA, boundary_events(kind, rowsB))
    stB = ps.n > 0 ? Int(ps.state[1]) : 0
    phB = ps.n > 0 ? ps.sph[1] - ph1 : NaN
    mB = isempty(rowsB) ? -1 : Int(rowsB[1].m)
    return (; evA, evT, stA, stB, phA, phB, mA, mB, rowsB)
end

# the photons of the emitter of ev_world over [w0, w1) with intensity 1 before te and hi from te on
switch_photons(te, hi, w0, w1) =
    te <= w0 ? 1000.0 * hi * (w1 - w0) : (te >= w1 ? 1000.0 * (w1 - w0) : 1000.0 * (te - w0) + 1000.0 * hi * (w1 - te))

@testset "closedloop/boundary_events_property" begin
    uex = UniformExcitation()
    wins = Dict(:interior => (0.004, 0.006, 0.008), :frame => (0.01, 0.02, 0.03))
    fwins = Dict(:interior => (0.004, 0.006, 0.008), :frame => (0.0, 0.01, 0.02))   # formation: a fresh world
    bdk = DimerKinetics(k_on=Inf, r_react=0.01, k_off=1.0, D_rot=0.0, d_dimer=0.005)
    ζ = randexp(copy(pair_world([bpop()], DimerKinetics(k_on=50.0, r_react=0.01, k_off=0.0, D_rot=0.0, d_dimer=0.005);
                                n_sub=1).rng))
    # recorded event: it takes effect in A exactly when te < a1, once, at te, inside the window that recorded it
    function check_recorded(kind, mode, win, te, r)
        a0, a1, b1 = win
        took = !isempty(r.evA)
        @test took == (te < a1)
        if took
            @test r.evA == [te]
            @test a0 <= te < a1
        end
        if hasproperty(r, :evT)
            @test length(r.evT) == 1
            rec = only(r.evT)
            (took || kind !== :bleach || te == a1) && @test rec == te   # a later bleach is recomputed in B from the leftover budget
            @test (took ? a0 <= rec < a1 : a1 <= rec < b1)
            @test 0.0 <= rec < (mode === :interior ? 0.01 : b1)
        end
    end
    for mode in (:interior, :frame)
        win = wins[mode]
        a0, a1, b1 = win
        targets = (prevfloat(a1), a1, nextfloat(a1))
        for te in targets
            @testset "departure $mode $te" begin
                r = boundary_run(() -> ev_world(), (w, ps) -> (ps.t_depart[1] = te; nothing), uex, :departure, mode, win)
                check_recorded(:departure, mode, win, te, r)
            end
            @testset "birth $mode $te" begin
                mk() = ev_world(lifetime=1e3, birth_rate=1e-9)
                r = boundary_run(mk, (w, ps) -> (ps.n = 0; ps.t_next_birth = te; nothing), uex, :birth, mode, win)
                check_recorded(:birth, mode, win, te, r)
            end
            @testset "switch $mode $te" begin
                exc = JumpExc(te, 1e18)
                r = boundary_run(() -> ev_world(), (w, ps) -> nothing, exc, :switch, mode, win)
                @test r.phA ≈ switch_photons(te, 1e18, a0, a1) rtol = 1e-9
                @test r.phB ≈ switch_photons(te, 1e18, a1, b1) rtol = 1e-9
                @test (r.phA > 100) == (te < a1)
            end
            @testset "break $mode $te" begin
                mk() = (w = pair_world([bpop()], bdk; n_sub=1); (w, w.pops[1]))
                function setup!(w, ps)
                    @test ps.partner[1] == 2 && ps.partner[2] == 1
                    ps.t_break_due[1] = te
                    ps.t_break_due[2] = te
                end
                r = boundary_run(mk, setup!, uex, :break, mode, win)
                check_recorded(:break, mode, win, te, r)
            end
            @testset "bleach $mode $te" begin
                f(b) = a0 + (0.0 + b / 1000.0)
                for (b, te′) in attain_all(f, 1000.0 * (te - a0), te)
                    r = boundary_run(() -> ev_world(budget=1e6), (w, ps) -> (ps.budget[1] = b; nothing), uex, :bleach, mode, win)
                    check_recorded(:bleach, mode, win, te′, r)
                end
            end
            @testset "exit $mode $te" begin
                f(c) = a0 + (0.0 + c / 1000.0)
                for (c, te′) in attain_all(f, 1000.0 * (te - a0), te)
                    mk() = ev_world(fluor=two_state(1000.0, 1000.0, 1e-9))
                    r = boundary_run(mk, (w, ps) -> (ps.state[1] = 1; ps.clock[1] = c; nothing), uex, :exit, mode, win)
                    @test r.stA == (te′ < a1 ? 2 : 1)
                    @test r.stB == 2
                end
            end
        end
        # formation: a fresh world, whose first draw of the sub-step is E
        fw = fwins[mode]
        for te in (prevfloat(fw[2]), fw[2], nextfloat(fw[2]))
            @testset "formation $mode $te" begin
                f(k) = fw[1] + ζ / k
                for (k, te′) in attain(f, ζ / (te - fw[1]), te; increasing=false)
                    dk = DimerKinetics(k_on=k, r_react=0.01, k_off=0.0, D_rot=0.0, d_dimer=0.005)
                    mk() = (w = pair_world([bpop()], dk; n_sub=1); (w, w.pops[1]))
                    r = boundary_run(mk, (w, ps) -> nothing, uex, :formation, mode, fw; warm=false, run_b=false)
                    took = !isempty(r.evA)
                    @test took == (te′ < fw[2])
                    took && @test r.evA == [te′]
                end
            end
        end
        # a bleach or state-1 exit due at a time <= a1 is never lost when the excitation is 0 from a1 on: before a1 it
        # fires in A, at exactly a1 it is left pending (budget or clock exactly 0) and fires at a1, after a1 it has a
        # positive residual and does not fire while dark
        @testset "pending at zero excitation $mode" begin
            exc = JumpExc(a1, 0.0)
            f(b) = a0 + (0.0 + b / 1000.0)
            for te in targets, (b, te′) in attain_all(f, 1000.0 * (te - a0), te)
                r = boundary_run(() -> ev_world(budget=1e6), (w, ps) -> (ps.budget[1] = b; nothing), exc, :pending_bleach, mode, win)
                if te′ < a1
                    @test r.evA == [te′] && r.evT == [te′]
                elseif te′ == a1
                    @test isempty(r.evA) && r.mA == 1
                    @test r.evT == [a1]
                    @test r.mB == 0
                else
                    @test isempty(r.evT) && r.mB == 1
                end
                mk() = ev_world(fluor=two_state(1000.0, 1000.0, 1e-9))
                r = boundary_run(mk, (w, ps) -> (ps.state[1] = 1; ps.clock[1] = b; nothing), exc, :exit, mode, win)
                @test r.stA == (te′ < a1 ? 2 : 1)
                @test r.stB == (te′ <= a1 ? 2 : 1)
            end
        end
        # a bleach after an earlier event of the same sub-step (a dark-to-lit exit, so τ > 0), beside a departure:
        # the earliest absolute time below a1 fires (ties: the departure), and a bleach due at a1 is pending
        @testset "bleach after an exit, beside a departure $mode" begin
            c0 = 1000.0 * 0.5 * (a1 - a0)
            f(b) = a0 + (c0 / 1000.0 + b / 1000.0)
            mk() = ev_world(fluor=two_state(1000.0, 1e-9, 1000.0), budget=1e6)
            for tt in targets, td in targets, (b, tb) in attain_all(f, 1000.0 * (tt - a0) - c0, tt)
                setup!(w, ps) = (ps.state[1] = 2; ps.clock[1] = c0; ps.budget[1] = b; ps.t_depart[1] = td; nothing)
                rd = boundary_run(mk, setup!, uex, :departure, mode, win)
                rb = boundary_run(mk, setup!, uex, :bleach, mode, win)
                ev = td <= tb ? (td < a1 ? :departure : :none) : (tb < a1 ? :bleach : :none)
                if ev === :departure
                    @test rd.evA == [td] && rd.evT == [td] && isempty(rb.evT)
                    @test a0 <= td < a1
                elseif ev === :bleach
                    @test rb.evA == [tb] && rb.evT == [tb] && isempty(rd.evT)
                    @test a0 <= tb < a1
                else
                    @test isempty(rd.evA) && isempty(rb.evA) && rd.stA == 1    # present, lit, at the end of A
                    if tb == a1                                   # pending: the first of the two at a1
                        td == a1 ? (@test rd.evT == [a1] && isempty(rb.evT)) : (@test rb.evT == [a1] && isempty(rd.evT))
                    end
                end
            end
        end
        # a state-1 exit after a declared switch, beside a bleach: the earliest below a1 fires (ties: the bleach),
        # and a pending one fires at a1
        @testset "exit after a switch, beside a bleach $mode" begin
            ts = a0 + 0.5 * (a1 - a0)
            exc = JumpExc(ts, 1.0)
            f(c) = a0 + ((ts - a0) + max((c - 1000.0 * (ts - a0)) / 1000.0, 0.0))
            mk() = ev_world(fluor=two_state(1000.0, 1000.0, 1e-9), budget=1e6)
            for tt in targets, tbt in targets, (c, tx) in attain_all(f, 1000.0 * (tt - a0), tt), (b, tb) in attain_all(f, 1000.0 * (tbt - a0), tbt)
                setup!(w, ps) = (ps.state[1] = 1; ps.clock[1] = c; ps.budget[1] = b; nothing)
                r = boundary_run(mk, setup!, exc, :bleach, mode, win)
                ev = min(tb, tx) < a1 ? (tb <= tx ? :bleach : :exit) : :none
                if ev === :bleach
                    @test r.evA == [tb] && r.evT == [tb] && r.stA == 0
                    @test a0 <= tb < a1
                elseif ev === :exit
                    @test isempty(r.evT) && r.stA == 2 && r.stB == 2    # dark after the exit: the bleach is lost
                else
                    @test isempty(r.evA) && r.stA == 1 && r.mA == 1
                    if tb == a1
                        @test r.evT == [a1]
                    elseif tx == a1
                        @test isempty(r.evT) && r.stB == 2
                    end
                end
            end
        end
    end
    # Codex's three cases
    @testset "three sub-steps, budget 10: the bleach at 0.01 is frame 2's" begin
        w, ps = ev_world(n_sub=3, budget=1e6)
        ps.budget[1] = 10.0
        SMLMSim.step!(w, 0.0, 0.01)
        r = only(frame_truth(w))
        @test r.m == 1 && isnan(r.t_bleach)
        SMLMSim.step!(w, 0.01, 0.02)
        r = only(frame_truth(w))
        @test r.m == 0
        @test r.t_bleach == 0.01
    end
    @testset "eight sub-steps over [0.03, 0.04), the budget of exactly that exposure" begin
        pop = Population(density=0.0, fluor=one_state(1000.0), psf=GaussianPSF(0.05), budget=1e6)
        w = SimWorld(StableRNG(1), cam32(), [pop]; n_sub=8, margin=0.0, t0=0.03)
        ps = w.pops[1]
        _add_emitter!(w, ps, 0.03)
        place!(ps, [1.6], [1.6])
        ps.budget[1] = 1000 * (0.04 - 0.03)
        SMLMSim.step!(w, 0.03, 0.04)
        r = only(frame_truth(w))
        @test r.m == 1 && isnan(r.t_bleach)
        SMLMSim.step!(w, 0.04, 0.05)
        r = only(frame_truth(w))
        @test r.m == 0
        @test 0.04 <= r.t_bleach < 0.05
    end
    @testset "a spot that switches off at 0.01 leaves the bleach pending" begin
        spot = Spot(x=1.6, y=1.6, σ=0.5, gain=1.0, z_R=Inf, t_off=0.01)
        sexc = SpotExcitation(base=0.0, spots=[spot])
        w, ps = ev_world(n_sub=1, budget=1e6)
        ps.budget[1] = 10.0
        SMLMSim.step!(w, 0.0, 0.01, sexc)
        @test only(frame_truth(w)).m == 1
        SMLMSim.step!(w, 0.01, 0.02, sexc)
        r = only(frame_truth(w))
        @test r.m == 0
        @test r.t_bleach == 0.01
    end
    # Codex's round-3 fixtures, one immobile emitter at one ulp of t
    ulpfl = GenericFluor(; γ=1.0, q=[-1e-9 1e-9; 1.0 -1.0])
    inside(w, rows) = all(r -> all(t -> isnan(t) || w.t_a <= t < w.t_b, (r.t_birth, r.t_bleach, r.t_depart, r.t_form, r.t_break)), rows)
    @testset "a rejected bleach at t1 does not mask a departure at prevfloat(t1)" begin
        t0 = 0.003
        t1 = t0 + 0.01
        w, ps = ev_world(n_sub=1, budget=1e6, fluor=ulpfl, t0=t0)
        ps.state[1] = 2
        ps.clock[1] = 0.006999999999999999
        ps.budget[1] = 0.0030000000000000005
        ps.t_depart[1] = prevfloat(t1)
        SMLMSim.step!(w, t0, t1)
        r = only(frame_truth(w))
        @test r.t_depart == prevfloat(t1)
        @test ps.n == 0
        @test inside(w, frame_truth(w))
        SMLMSim.step!(w, t1, t1 + 0.01)
        @test isempty(frame_truth(w))
    end
    @testset "a bleach at a time below t1 is not deferred for a relative delay above the remainder" begin
        (t0, t1, c, b) = (0.09738069177455938, 0.10192574722038701, 0.0041244712588255, 0.0004205841870021238)
        w, ps = ev_world(n_sub=1, budget=1e6, fluor=ulpfl, t0=t0)
        ps.state[1] = 2
        ps.clock[1] = c
        ps.budget[1] = b
        SMLMSim.step!(w, t0, t1)
        r = only(frame_truth(w))
        @test t0 + (c + b) < t1
        @test r.t_bleach == t0 + (c + b)
        @test r.m == 0
        @test inside(w, frame_truth(w))
    end
    @testset "a bleach or exit due at t1 survives a switch-off at t1" begin
        spot = Spot(x=1.6, y=1.6, σ=0.5, gain=1.0, z_R=Inf, t_off=0.03)
        sexc = SpotExcitation(base=0.0, spots=[spot])
        w, ps = ev_world(n_sub=1, budget=1e6, t0=0.02)         # one-state, budget 10 at 1000/s: due at 0.02 + 0.01
        ps.budget[1] = 10.0
        SMLMSim.step!(w, 0.02, 0.03, sexc)
        @test only(frame_truth(w)).m == 1
        SMLMSim.step!(w, 0.03, 0.04, sexc)
        r = only(frame_truth(w))
        @test r.m == 0
        @test r.t_bleach == 0.03
        @test inside(w, frame_truth(w))
        w, ps = ev_world(n_sub=1, fluor=GenericFluor(; γ=1000.0, q=[-100.0 100.0; 1e-12 -1e-12]), t0=0.02)
        ps.state[1] = 1
        ps.clock[1] = 1.0
        SMLMSim.step!(w, 0.02, 0.03, sexc)
        @test only(frame_truth(w)).lit == 1 && ps.state[1] == 1
        SMLMSim.step!(w, 0.03, 0.04, sexc)
        @test ps.state[1] == 2
        @test only(frame_truth(w)).lit == 0
    end
end
