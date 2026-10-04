using SMLMSim, Test, Statistics, Random, StableRNGs, MicroscopePSFs
using SMLMSim: frame_truth
using SMLMSim.Stepper: _add_emitter!

cam32() = IdealCamera(1:32, 1:32, 0.1)
one_state(γ) = GenericFluor(; γ, q=zeros(1, 1))
two_state(γ, k_off, k_on) = GenericFluor(; γ, q=[-k_off k_off; k_on -k_on])

# largest distance between the empirical CDF of ts and the exponential CDF of the given mean
function ks_exp(ts, mean)
    s = sort(ts)
    n = length(s)
    return maximum(max(abs(k / n - (1 - exp(-s[k] / mean))), abs((k - 1) / n - (1 - exp(-s[k] / mean))))
                   for k in 1:n)
end

# the plain Gaussian-beam law, written inline in spot order
function plain_law(base, spots, x, y, z)
    I = base
    for (xk, yk, σ, gain, z_R) in spots
        s = σ * sqrt(1 + (z / z_R)^2)
        r2 = (x - xk)^2 + (y - yk)^2
        I += gain * (σ / s)^2 * exp(-r2 / (2 * s^2))
    end
    return I
end

@testset "closedloop/spot_tilt_evan" begin
    σ, gain, base = 0.0934, 5.0, 1.0
    # tilt: the maximum at z = 0.5 sits at (x_k, y_k) + z tan θ (cos φ, sin φ)
    for φ in (0.0, 0.7, 2.5)
        e = SpotExcitation(; base, spots=[Spot(; x=1.6, y=1.7, σ, gain, z_R=0.23)], tilt=(2.0, φ))
        cx, cy = 1.6 + 1.0 * cos(φ), 1.7 + 1.0 * sin(φ)
        grid = (-300:300) .* 1e-3
        best, bx, by = -Inf, 0.0, 0.0
        for dx in grid, dy in grid
            I = e(cx + dx, cy + dy, 0.5, 0.0)
            I > best && ((best, bx, by) = (I, cx + dx, cy + dy))
        end
        @test hypot(bx - cx, by - cy) <= 1e-3
    end
    # evanescent split: gain f at the glass, 1 - f in the beam; each term checked apart
    d_evan = 0.1
    spots = [Spot(; x=1.6, y=1.6, σ, gain, z_R=Inf)]
    e0 = SpotExcitation(; base, spots)
    e2 = SpotExcitation(; base, spots, f_evan=0.2, d_evan)
    ev(z) = e2(1.6, 1.6, z, 0.0) - base - 0.8 * (e0(1.6, 1.6, z, 0.0) - base)
    @test isapprox(e0(1.6, 1.6, 0.0, 0.0) - base, gain; rtol=1e-12)
    @test isapprox(e2(1.6, 1.6, 0.0, 0.0) - base - ev(0.0), 0.8 * gain; rtol=1e-12)
    @test isapprox(ev(0.0), 0.2 * gain; rtol=1e-12)
    @test ev(3 * d_evan) < 0.05 * ev(0.0)
    @test isapprox(ev(d_evan), 0.2 * gain * exp(-1); rtol=1e-12)
    @test e2(1.6, 1.6, -0.1, 0.0) - base ≈ 0.8 * (e0(1.6, 1.6, -0.1, 0.0) - base)   # no evanescent term below the glass
    # defaults reproduce the plain law bitwise
    ps = [(1.0, 1.2, 0.09, 4.0, 0.23), (1.9, 1.5, 0.12, 30.0, 0.5), (1.5, 2.2, 0.2, 2.0, Inf)]
    e = SpotExcitation(; base=0.7, spots=[Spot(; x, y, σ, gain, z_R) for (x, y, σ, gain, z_R) in ps])
    r = MersenneTwister(5)
    same = true
    for _ in 1:2000
        x, y, z = 3 * rand(r), 3 * rand(r), 2 * rand(r) - 0.5
        same &= e(x, y, z, 0.0) === plain_law(0.7, ps, x, y, z)
    end
    @test same
end

@testset "closedloop/spot_keywords" begin
    ok = (x=1.0, y=1.0, σ=0.0934, gain=5.0, z_R=0.23)
    @test Spot(; ok...) isa Spot
    for k in keys(ok)
        @test_throws UndefKeywordError Spot(; (; (j => ok[j] for j in keys(ok) if j != k)...)...)
    end
    @test_throws MethodError Spot(1.0, 1.0, 0.0934, 5.0, 0.23)
    @test_throws ArgumentError Spot(; ok..., σ=0.0)
    @test_throws ArgumentError Spot(; ok..., gain=-1.0)
    @test_throws ArgumentError Spot(; ok..., z_R=0.0)
    @test Spot(; ok..., z_R=Inf).z_R == Inf
    @test Spot(; ok...).t_on == -Inf && Spot(; ok...).t_off == Inf
    @test_throws ArgumentError Spot(; ok..., t_on=1.0, t_off=1.0)
    @test_throws ArgumentError Spot(; ok..., t_on=2.0, t_off=1.0)
    @test_throws ArgumentError Spot(; ok..., t_on=NaN)
    @test_throws ArgumentError Spot(; ok..., t_off=NaN)
    @test Spot(; ok..., t_on=-Inf, t_off=Inf) isa Spot
    # the stored Float64 values are checked, not the inputs: a value that passes before conversion cannot store
    # a broken one
    @test_throws ArgumentError Spot(; ok..., σ=big"1e-400")
    @test_throws ArgumentError Spot(; ok..., t_on=big(2)^53, t_off=big(2)^53 + 1)
    @test_throws ArgumentError Spot(; ok..., σ=Inf)
    @test_throws ArgumentError Spot(; ok..., gain=Inf)
    @test_throws ArgumentError Spot(; ok..., x=NaN)
    @test_throws ArgumentError Spot(; ok..., y=Inf)
    sp = [Spot(; ok...)]
    @test_throws UndefKeywordError SpotExcitation()
    @test_throws MethodError SpotExcitation(1.0, sp)
    @test_throws ArgumentError SpotExcitation(; spots=sp, base=-0.1)
    @test_throws ArgumentError SpotExcitation(; spots=sp, f_evan=-0.1)
    @test_throws ArgumentError SpotExcitation(; spots=sp, f_evan=1.1)
    @test_throws ArgumentError SpotExcitation(; spots=sp, d_evan=0.0)
    @test_throws ArgumentError SpotExcitation(; spots=sp, λ=0.0)
    @test_throws ArgumentError SpotExcitation(; spots=sp, n=0.0)
    @test_throws ArgumentError SpotExcitation(; spots=sp, tilt=(-1.0, 0.0))
    @test_throws ArgumentError SpotExcitation(; spots=sp, tilt=(Inf, 0.0))
    @test_throws ArgumentError SpotExcitation(; spots=sp, tilt=(NaN, 0.0))
    @test_throws ArgumentError SpotExcitation(; spots=sp, tilt=(big"1e400", 0.0))
    @test_throws ArgumentError SpotExcitation(; spots=sp, base=Inf)
    @test_throws ArgumentError SpotExcitation(; spots=sp, λ=Inf)
    @test_throws ArgumentError SpotExcitation(; spots=sp, n=NaN)
    @test SpotExcitation(; spots=sp, base=0.0, f_evan=1.0) isa SpotExcitation
    # the rig spot does not warn; a z_R 7x too long does, once per construction, naming the spot
    @test_logs SpotExcitation(; spots=sp)
    far = [Spot(; ok..., z_R=1.7)]
    @test_logs (:warn, r"spot 1: z_R = 1.7") SpotExcitation(; spots=far)
    @test_logs (:warn, r"spot 1: z_R = 1.7") SpotExcitation(; spots=far)
    @test_logs (:warn, r"spot 2: z_R = 0.01") SpotExcitation(; spots=[Spot(; ok...), Spot(; ok..., z_R=0.01)])
    @test_logs SpotExcitation(; spots=[Spot(; ok..., z_R=Inf)])
end

@testset "closedloop/bleach_time_spot" begin
    B, γ, T = 2000.0, 1e5, 1.0
    σ, gain, base, z_R = 0.0934, 4.0, 1.0, 0.23
    exc = SpotExcitation(; base, spots=[Spot(; x=1.6, y=1.6, σ, gain, z_R)])
    for (z, I0) in ((0.0, base + gain), (z_R, base + gain / 2))
        pop = Population(density=0.0, fluor=one_state(γ), budget=B, z=(z, z), psf=GaussianPSF(0.05))
        w = SimWorld(StableRNG(8), cam32(), [pop]; n_sub=1, margin=0.0)
        ps = w.pops[1]
        for _ in 1:5000
            _add_emitter!(w, ps, 0.0)
        end
        fill!(view(ps.x, 1:5000), 1.6)
        fill!(view(ps.y, 1:5000), 1.6)
        SMLMSim.step!(w, 0.0, T, exc)
        tb = [r.t_bleach for r in frame_truth(w) if !isnan(r.t_bleach)]
        @test length(tb) == 5000
        @test ks_exp(tb, B / (γ * I0)) * sqrt(length(tb)) < 1.63   # p > 0.01
    end
end

# one immobile, one-state emitter at the spot centre over the frames 1..nframes of T = 10 ms; returns the
# last frame's truth row and the sum of the signal map
function one_emitter_frame(spot, base, nframes; γ=1e5, T=0.01)
    pop = Population(density=0.0, fluor=one_state(γ), budget=Inf, psf=GaussianPSF(0.05))
    w = SimWorld(StableRNG(9), cam32(), [pop]; n_sub=8, margin=0.0)
    _add_emitter!(w, w.pops[1], 0.0)
    w.pops[1].x[1] = 1.6
    w.pops[1].y[1] = 1.6
    exc = SpotExcitation(; base, spots=[spot])
    for k in 1:nframes
        SMLMSim.step!(w, (k - 1) * T, k * T, exc)
    end
    return only(frame_truth(w)), sum(w.signal)
end

# The expectation recorded by the stepper, not a noisy draw: FrameTruth's photons (the noise-free emitted
# total) and the signal map's sum. `lit` is the fraction of T in the emitting state, which excitation does
# not change; the switch shows in `excitation`, the presence-weighted mean intensity. `lit` is `FrameTruth`'s
# emitting-state fraction (state 1 with `m > 0`, see its docstring), so it is 1 here whatever the spot does.
@testset "closedloop/spot_switch_mid_exposure" begin
    γ, T, σ, gain, z_R = 1e5, 0.01, 0.0934, 3.0, 0.23
    for nframes in (1, 3), base in (0.0, 1.0)
        t_a, t_b = (nframes - 1) * T, nframes * T
        # on at 0.37 of the way into the exposure: not a sub-step boundary (n_sub = 8)
        t_on = t_a + 0.37 * T
        row, tot = one_emitter_frame(Spot(; x=1.6, y=1.6, σ, gain, z_R, t_on), base, nframes)
        lit_time = t_b - t_on
        @test isapprox(row.photons, γ * (base * T + gain * lit_time); rtol=1e-9)
        @test isapprox(tot, row.photons; rtol=2e-6)
        @test isapprox(row.excitation, base + gain * lit_time / T; rtol=1e-9)
        @test isapprox(row.lit, 1; rtol=1e-12)
        # off at 0.62 of the way in
        t_off = t_a + 0.62 * T
        row, tot = one_emitter_frame(Spot(; x=1.6, y=1.6, σ, gain, z_R, t_off), base, nframes)
        @test isapprox(row.photons, γ * (base * T + gain * (t_off - t_a)); rtol=1e-9)
        @test isapprox(tot, row.photons; rtol=2e-6)
        @test isapprox(row.excitation, base + gain * (t_off - t_a) / T; rtol=1e-9)
        # a window inside the exposure
        row, _ = one_emitter_frame(Spot(; x=1.6, y=1.6, σ, gain, z_R, t_on=t_a + 0.21 * T, t_off=t_a + 0.84 * T), base, nframes)
        @test isapprox(row.photons, γ * (base * T + gain * 0.63 * T); rtol=1e-9)
    end
    # a spot that is off for the whole exposure adds nothing, and one that is always on adds gain
    row, _ = one_emitter_frame(Spot(; x=1.6, y=1.6, σ, gain, z_R, t_on=5.0), 1.0, 2)
    @test isapprox(row.photons, γ * T; rtol=1e-9)
    row, _ = one_emitter_frame(Spot(; x=1.6, y=1.6, σ, gain, z_R), 0.0, 2)
    @test isapprox(row.photons, γ * gain * T; rtol=1e-9)
    # right-continuous: on from t_on, off from t_off
    e = SpotExcitation(; base=0.0, spots=[Spot(; x=0.0, y=0.0, σ, gain, z_R, t_on=0.1, t_off=0.2)])
    @test e(0.0, 0.0, 0.0, prevfloat(0.1)) == 0.0
    @test e(0.0, 0.0, 0.0, 0.1) == gain
    @test e(0.0, 0.0, 0.0, prevfloat(0.2)) == gain
    @test e(0.0, 0.0, 0.0, 0.2) == 0.0
end

# one immobile emitter at the centre of a spot that turns on at 3.7 ms of a 10 ms exposure (I = 0 before, 3 after),
# with a finite budget or a state-1 exit that the switch must move
@testset "closedloop/spot_switch_budget_clock" begin
    γ, T, t_on = 1e5, 0.01, 0.0037
    exc = SpotExcitation(; base=0.0, spots=[Spot(; x=1.6, y=1.6, σ=0.0934, gain=3.0, z_R=0.23, t_on)])
    function frame(fluor, prep!)
        pop = Population(density=0.0, fluor=fluor, budget=1e9, psf=GaussianPSF(0.05))
        w = SimWorld(StableRNG(11), cam32(), [pop]; n_sub=8, margin=0.0)
        ps = w.pops[1]
        _add_emitter!(w, ps, 0.0)
        ps.x[1] = 1.6
        ps.y[1] = 1.6
        prep!(ps)
        SMLMSim.step!(w, 0.0, T, exc)
        return only(frame_truth(w))
    end
    # (a) the budget runs out 300/(γ·3) after the switch
    row = frame(one_state(γ), ps -> (ps.budget[1] = 300.0))
    @test isapprox(row.t_bleach, t_on + 300 / (γ * 3); rtol=1e-9)
    @test isapprox(row.photons, 300; rtol=1e-9)
    @test row.m == 0
    # (b) the state-1 exit clock runs only while lit: it leaves 1/(100·3) s after the switch
    row = frame(GenericFluor(; γ, q=[-100.0 100.0; 1e-12 -1e-12]), ps -> (ps.state[1] = 1; ps.clock[1] = 1.0))
    @test isapprox(row.photons, 1000; rtol=1e-9)
    @test isapprox(row.lit, (t_on + 1 / 300) / T; rtol=1e-9)
end

@testset "closedloop/spot_next_switch" begin
    mk(; kw...) = Spot(; x=0.0, y=0.0, σ=0.1, gain=1.0, z_R=0.23, kw...)
    e = SpotExcitation(; spots=[mk(t_on=0.1, t_off=0.3), mk(t_on=0.2, t_off=0.3), mk(), mk(t_on=0.5)])
    @test next_switch(e, -Inf) == 0.1
    @test next_switch(e, -1.0) == 0.1
    @test next_switch(e, 0.0) == 0.1
    @test next_switch(e, 0.1) == 0.2     # at a switch: the next one after it
    @test next_switch(e, prevfloat(0.1)) == 0.1
    @test next_switch(e, 0.15) == 0.2
    @test next_switch(e, 0.2) == 0.3
    @test next_switch(e, 0.3) == 0.5     # a tie at 0.3 counts once
    @test next_switch(e, 0.5) == Inf
    @test next_switch(e, 7.0) == Inf
    @test next_switch(e, Inf) == Inf
    @test next_switch(SpotExcitation(; spots=[mk()]), 0.0) == Inf
    @test next_switch(SpotExcitation(; spots=Spot[]), 0.0) == Inf
    @test next_switch(SpotExcitation(; spots=[mk(t_off=0.7)]), 0.7) == Inf
    @test next_switch(SpotExcitation(; spots=[mk(t_off=0.7)]), 0.6) == 0.7
    @test next_switch(SpotExcitation(; spots=[mk(t_on=-Inf, t_off=Inf)]), -Inf) == Inf
    # the allocation-free loop (measured inside a function: at top level of a testset Julia 1.10 boxes the result)
    ns_alloc(e, t) = @allocated next_switch(e, t)
    ns_alloc(e, 0.1)
    @test ns_alloc(e, 0.1) == 0
end

@testset "closedloop/spot_zero_alloc" begin
    pop = Population(density=3.0, fluor=two_state(1e4, 300.0, 200.0), budget=3000.0, lifetime=0.5,
                     mobility=[(1.0, 1.0)], z=(-0.1, 0.3), psf=GaussianPSF(0.05))
    w = SimWorld(StableRNG(23), cam32(), [pop]; n_sub=4, margin=0.0)
    mk(x, y; kw...) = Spot(; x, y, σ=0.0934, gain=20.0, z_R=0.23, kw...)
    # steady spots, and spots that switch inside the measured exposures (t from 2.0 to 3.0 s)
    exc = SpotExcitation(; base=1.0, spots=[mk(0.8, 0.8), mk(2.0, 2.4; t_on=2.013, t_off=2.467),
                                            mk(2.6, 0.9; t_on=2.801)], f_evan=0.2, tilt=(0.5, 1.0))
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
