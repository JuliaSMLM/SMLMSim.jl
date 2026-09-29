using SMLMSim, Test, Statistics, Random, StableRNGs, MicroscopePSFs
using SMLMSim: frame_truth, DimerKinetics, DiffusionSMLMConfig, FrameTruth
using SMLMSim.Stepper: _add_emitter!
using Random: randexp

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

const DIMER_GOLDEN = 55106.829324442115

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
    function vis_world(xs, ys; pops=nothing, merge_radius=0.0, mult=(1, 1))
        pops === nothing && (pops = [Population(name=:A, density=0.0, fluor=fl, psf=GaussianPSF(0.05)),
                                     Population(name=:B, density=0.0, fluor=fl, multiplicity=mult[2], psf=GaussianPSF(0.05))])
        w = SimWorld(StableRNG(7), cam32(), pops; n_sub=4, margin=0.0, dimers=dk, merge_radius)
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
end
