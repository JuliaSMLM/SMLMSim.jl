using SMLMSim, Test, Statistics, Random, StableRNGs, MicroscopePSFs
using SMLMSim: frame_truth
using SMLMSim.Stepper: _add_emitter!

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
            for t in (r.t_birth, r.t_depart, r.t_bleach)
                @test isnan(t) || (ta <= t < tb)
            end
            (isnan(r.t_depart) && isnan(r.t_bleach)) || push!(ended, r.id)
            @test 0 <= r.lit <= 1
            @test r.frame == k
        end
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
