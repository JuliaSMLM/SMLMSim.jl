using SMLMSim, Test, Random, Statistics

# Directed motion: a constant velocity with diffusion on top, for every monomer, per population
# (monomer_velocity) or per molecule (track_velocity, track_D).

const DM = SMLMSim.InteractionDiffusion

# Every step recorded (exposure = frame period = dt), no binding (r_react tiny, low density)
dm_config(; kw...) = DiffusionSMLMConfig(; density=1.0, box_size=10.0, diff_monomer=0.05, diff_dimer=0.0,
    r_react=1e-6, dt=0.01, t_max=0.5, camera_framerate=100.0, camera_exposure=0.01, kw...)

coords(e) = hasproperty(e, :z) ? (e.x, e.y, e.z) : (e.x, e.y)

# Per-step displacements of monomer steps (both ends :monomer, one dt apart), by track, at the minimum image
function monomer_steps(smld, box; dt=0.01)
    by = Dict{Int,Vector{Any}}()
    for e in smld.emitters
        push!(get!(by, e.track_id, Any[]), e)
    end
    steps = Dict{Int,Vector{Vector{Float64}}}()
    for (id, r) in by
        sort!(r, by=e -> e.timestamp)
        for i in 1:length(r)-1
            a, b = r[i], r[i+1]
            (a.state == :monomer && b.state == :monomer && isapprox(b.timestamp - a.timestamp, dt; rtol=1e-6)) || continue
            d = [Float64(q) - Float64(p) for (p, q) in zip(coords(a), coords(b))]
            push!(get!(steps, id, Vector{Float64}[]), d .- box .* round.(d ./ box))
        end
    end
    return steps
end

# Mean within 5 standard errors of v*dt and variance within 5 standard errors of 2*D*dt, on every axis
function check_stats(ds, v, D; dt=0.01)
    @test length(ds) > 1000
    n = length(ds)
    for k in eachindex(v)
        x = [d[k] for d in ds]
        σ2 = 2 * D * dt
        @test abs(mean(x) - v[k] * dt) < 5 * sqrt(σ2 / n)
        @test abs(var(x) - σ2) < 5 * σ2 * sqrt(2 / (n - 1))
    end
end

@testset "Directed motion" begin
    @testset "one velocity for every monomer, 2D and 3D" begin
        for (nd, v) in ((2, (0.8, -0.3)), (3, (0.4, 0.2, -0.6)))
            Random.seed!(9100 + nd)
            p = dm_config(ndims=nd, monomer_velocity=[v])
            smld, _ = simulate(p; override_count=400, γ=1000.0)
            steps = monomer_steps(smld, p.box_size)
            check_stats(reduce(vcat, values(steps)), v, p.diff_monomer)
            # a run with a velocity saves every track's
            mv = smld.metadata["monomer_velocity"]
            @test length(mv) == 400 && all(==(v), values(mv))
        end
    end

    @testset "pure directed motion: exact steps, no diffusion" begin
        Random.seed!(9110)
        v = (0.5, 0.25, -1.0)
        p = dm_config(ndims=3, diff_monomer=0.0, monomer_velocity=[collect(v)])  # a vector entry is accepted
        @test p.monomer_velocity == [v]
        smld, _ = simulate(p; override_count=50, γ=1000.0)
        for d in reduce(vcat, values(monomer_steps(smld, p.box_size)))
            @test all(isapprox.(d, v .* p.dt; atol=1e-12))
        end
    end

    @testset "per population" begin
        Random.seed!(9120)
        mix = [(0.3, 0.0), (0.3, 0.02), (0.4, 0.1)]
        vel = [(0.0, 0.0), (1.0, 0.5), (-0.4, 0.0)]
        p = dm_config(monomer_mobility=mix, monomer_velocity=vel)
        smld, _ = simulate(p; override_count=900, γ=1000.0)
        cls, mv = smld.metadata["monomer_class"], smld.metadata["monomer_velocity"]
        @test length(mv) == 900 && all(mv[id] == vel[cls[id]] for id in keys(cls))
        steps = monomer_steps(smld, p.box_size)
        for k in 1:3
            ds = reduce(vcat, [s for (id, s) in steps if cls[id] == k])
            if k == 1
                @test all(d -> all(iszero, d), ds)  # immobile population: no drift, no diffusion
            else
                check_stats(ds, vel[k], mix[k][2])
            end
        end
    end

    @testset "per molecule: track_velocity and track_D" begin
        Random.seed!(9130)
        p = dm_config()
        tv = Dict(id => (1.0, -0.5) for id in 1:200)
        tD = Dict(id => 0.0 for id in 1:100)
        smld, _ = simulate(p; override_count=300, γ=1000.0, track_velocity=tv, track_D=tD)
        steps = monomer_steps(smld, p.box_size)
        directed = reduce(vcat, [steps[id] for id in 1:100])
        @test all(d -> all(isapprox.(d, (0.01, -0.005); atol=1e-12)), directed)
        check_stats(reduce(vcat, [steps[id] for id in 101:200]), (1.0, -0.5), p.diff_monomer)
        check_stats(reduce(vcat, [steps[id] for id in 201:300]), (0.0, 0.0), p.diff_monomer)
        mv = smld.metadata["monomer_velocity"]
        @test sort!(collect(keys(mv))) == 1:300 && all(mv[id] == (0.0, 0.0) for id in 201:300)
        @test sort!(collect(keys(smld.metadata["monomer_D"]))) == 1:100
        # the keywords leave the other molecules' population draws unchanged
        q = dm_config(monomer_mobility=[(0.5, 0.01), (0.5, 0.1)])
        Random.seed!(9131); a, _ = simulate(q; override_count=50, γ=1000.0)
        Random.seed!(9131); b, _ = simulate(q; override_count=50, γ=1000.0, track_D=Dict(3 => 0.5),
                                            track_velocity=Dict(4 => (1.0, 0.0)))
        @test a.metadata["monomer_class"] == b.metadata["monomer_class"]
        @test b.metadata["monomer_D"][3] == 0.5 && b.metadata["monomer_velocity"][4] == (1.0, 0.0)
    end

    @testset "no velocity: unchanged output" begin
        for kw in ((monomer_velocity=[(0.0, 0.0)],), (monomer_mobility=[(0.5, 0.0), (0.5, 0.1)],
                                                      monomer_velocity=[(0.0, 0.0), (0.0, 0.0)]))
            base = (; density=20.0, box_size=1.0, t_max=0.3, r_react=0.03, d_dimer=0.04, k_off=5.0)
            mix = haskey(kw, :monomer_mobility) ? (; monomer_mobility=kw.monomer_mobility) : (;)
            Random.seed!(9140); a, _ = simulate(DiffusionSMLMConfig(; base..., mix...); γ=1000.0)
            Random.seed!(9140); b, _ = simulate(DiffusionSMLMConfig(; base..., kw...); γ=1000.0)
            @test any(e -> e.state == :dimer, a.emitters)
            @test a.emitters == b.emitters
            @test !haskey(a.metadata, "monomer_velocity") && haskey(b.metadata, "monomer_velocity")
        end
    end

    @testset "no velocity: main's output and random stream (golden values)" begin
        # Recorded by running these three cases, as written (with this file's coords and DM), on main b016876: a
        # `git archive` of it developed in a temporary environment that adds Random. Identical on Julia 1.10, 1.11,
        # 1.12 and 1.13. Julia does not promise its random streams across versions, and a few rare randn draws differ
        # between them, so these cases are ones that agree on all four. The next rand() after each run pins the
        # random stream: an extra draw fails it. Re-record only when a toolchain change moves the random numbers.
        fingerprint(s) = (n = length(s.emitters);
            (n, count(e -> e.state == :dimer, s.emitters), [coords(s.emitters[i]) for i in (1, n ÷ 2, n)],
             sum(Float64(e.x) for e in s.emitters), rand()))
        Random.seed!(4711)
        s, _ = simulate(DiffusionSMLMConfig(density=20.0, box_size=1.0, t_max=0.3, r_react=0.03, d_dimer=0.04, k_off=5.0,
                                            monomer_mobility=[(0.5, 0.0), (0.5, 0.1)]); γ=1000.0)
        @test fingerprint(s) == (600, 178, [(0.3786698976711056, 0.23892995154241714),
                                            (0.09777886727368046, 0.28861049568788477),
                                            (0.09777886727368046, 0.28861049568788477)],
                                 264.1249680481553, 0.14501070246701409)
        Random.seed!(4713)
        s, _ = simulate(DiffusionSMLMConfig(ndims=3, boundary="reflecting", density=5.0, box_size=1.0, t_max=0.2,
                                            r_react=0.05, d_dimer=0.06, k_off=3.0, pair_mobility=:min,
                                            monomer_mobility=[(0.3, 0.0), (0.7, 0.2)]); γ=1000.0)
        @test fingerprint(s) == (100, 0, [(0.2832631894221316, 0.704742513176353, 0.4820392062953146),
                                          (0.3929977534116499, 0.5202338658078117, 0.9792805052040606),
                                          (0.5834611575571461, 0.5314309101411046, 0.8487075602517744)],
                                 48.66090921247954, 0.4164589717507211)
        Random.seed!(4713)
        start = [DM.DiffusingEmitter2D{Float32}(Float32(0.1 * i), Float32(0.05 * i), 10.0, 0.0, 1, 1, i, :monomer, nothing)
                 for i in 1:9]
        s, _ = simulate(DiffusionSMLMConfig(boundary="reflecting", box_size=1.0, t_max=0.3, r_react=0.04, d_dimer=0.05,
                                            diff_monomer=0.2, k_off=4.0); starting_conditions=start, γ=1000.0)
        @test fingerprint(s) == (270, 88, [(0.1f0, 0.05f0), (0.96876913f0, 0.75858253f0), (0.8480124f0, 0.7458643f0)],
                                 145.99921490997076, 0.8814223011579213)
    end

    @testset "bound pair" begin
        pair(; d=0.05) = [DM.DiffusingEmitter2D{Float64}(1.0, 1.0, 10.0, 0.0, 1, 1, 1, :dimer, 2),
                          DM.DiffusingEmitter2D{Float64}(1.0 + d, 1.0, 10.0, 0.0, 1, 1, 2, :dimer, 1)]
        run(pm; kw...) = (Random.seed!(9150); simulate(dm_config(diff_monomer=0.1, diff_dimer=0.0, diff_dimer_rot=0.0,
            k_off=0.0, d_dimer=0.05, pair_mobility=pm, t_max=0.2); starting_conditions=pair(), kw...)[1])
        track(s, id) = sort([e for e in s.emitters if e.track_id == id], by=e -> e.timestamp)
        disp(s, id) = (r = track(s, id); [(r[i+1].x - r[i].x, r[i+1].y - r[i].y) for i in 1:length(r)-1])
        # :fixed, no pair diffusion: the pair moves at the mean of its partners' velocities
        s = run(:fixed; track_velocity=Dict(1 => (1.0, 0.0), 2 => (0.0, 0.4)))
        @test all(e -> e.state == :dimer, s.emitters)
        for id in 1:2
            @test all(d -> isapprox(d[1], 0.005; atol=1e-12) && isapprox(d[2], 0.002; atol=1e-12), disp(s, id))
        end
        @test all(isapprox(hypot(a.x - b.x, a.y - b.y), 0.05; atol=1e-12) for (a, b) in zip(track(s, 1), track(s, 2)))
        # :min, a directed member with D = 0 is not immobile: the pair (D_pair = 0) only drifts
        s = run(:min; track_velocity=Dict(1 => (1.0, 0.0)), track_D=Dict(1 => 0.0))
        @test all(d -> isapprox(d[1], 0.005; atol=1e-12) && d[2] == 0, disp(s, 2))
        @test !any(r -> r.mixed, frame_dimer_truth(s))
        # :min, an immobile member (D = 0, no velocity) anchors the pair: its partner's velocity does not move it
        s = run(:min; track_velocity=Dict(1 => (1.0, 0.0)), track_D=Dict(2 => 0.0))
        @test all(d -> d == (0.0, 0.0), disp(s, 1)) && all(d -> d == (0.0, 0.0), disp(s, 2))
        @test all(r -> r.mixed, frame_dimer_truth(s))
        # :fixed warns when it moves an immobile member, and not when the D = 0 member drifts
        @test_logs (:warn, r"moved immobile") match_mode=:any run(:fixed; track_D=Dict(2 => 0.0),
                                                                  track_velocity=Dict(1 => (1.0, 0.0)))
        @test_logs run(:fixed; track_D=Dict(1 => 0.0, 2 => 0.0), track_velocity=Dict(1 => (1.0, 0.0), 2 => (1.0, 0.0)))
    end

    @testset "r_react = 0: no binding" begin
        # four D 0 movers through one point at the same moment, a crossing scene of a tracking benchmark: (4, 4)
        # at t = 1.5 s
        paths = [(1.0, 1.0, 2.0, 2.0), (7.0, 1.0, -2.0, 2.0), (1.0, 7.0, 2.0, -2.0), (7.0, 7.0, -2.0, -2.0)]
        start = [DM.DiffusingEmitter2D{Float64}(x, y, 10.0, 0.0, 1, 1, id, :monomer, nothing)
                 for (id, (x, y, _, _)) in enumerate(paths)]
        tv = Dict(id => (vx, vy) for (id, (_, _, vx, vy)) in enumerate(paths))
        cross(r) = (Random.seed!(9190); simulate(dm_config(r_react=r, t_max=3.0); starting_conditions=start, γ=1000.0,
                                                 track_velocity=tv, track_D=Dict(id => 0.0 for id in 1:4))[1])
        s = cross(0.0)
        @test all(e -> e.state == :monomer, s.emitters)
        # each passes straight through: start + v·t at every record
        for e in s.emitters
            x0, y0, vx, vy = paths[e.track_id]
            @test isapprox(e.x, x0 + vx * e.timestamp; atol=1e-9) && isapprox(e.y, y0 + vy * e.timestamp; atol=1e-9)
        end
        mid = [e for e in s.emitters if isapprox(e.timestamp, 1.5; atol=1e-9)]
        @test length(mid) == 4 && all(e -> isapprox(e.x, 4.0; atol=1e-9) && isapprox(e.y, 4.0; atol=1e-9), mid)
        # the same paths with r_react > 0 bind
        @test any(e -> e.state == :dimer, cross(0.01).emitters)
        @test DiffusionSMLMConfig(r_react=0.0).r_react == 0.0
        @test_throws ArgumentError DiffusionSMLMConfig(r_react=-0.01)
    end

    @testset "reflecting wall" begin
        Random.seed!(9160)
        p = dm_config(diff_monomer=0.0, box_size=2.0, boundary="reflecting", monomer_velocity=[(2.0, 0.0)], t_max=1.5)
        smld, _ = simulate(p; override_count=20, γ=1000.0)
        @test all(e -> 0 <= e.x <= p.box_size, smld.emitters)
        last_x = [e.x for e in smld.emitters if e.frame == smld.n_frames]
        @test all(x -> x >= p.box_size - 2.0 * p.dt - 1e-12, last_x)  # held within one step of the wall

        # A step longer than the box folds back off the walls as often as it crosses them (Codex review r2: one
        # reflection then a clamp left such a step at the wall), with drift and boundary in Float64
        mono(x; T=Float64) = DM.DiffusingEmitter2D{T}(T(x), T(0.5), 10.0, 0.0, 1, 1, 1, :monomer, nothing)
        box1(b) = dm_config(diff_monomer=0.0, box_size=1.0, boundary=b, r_react=0.0)
        step1(b, x, vx; T=Float64) = only(DM.update_system([mono(x; T)], box1(b), 0.01; track_v=Dict(1 => (vx, 0.0))))
        for T in (Float64, Float32)
            @test isapprox(step1("reflecting", 0.25, 300.0; T).x, 0.75; atol=1e-6)   # 3.25: off 1, off 0, off 1
            @test isapprox(step1("reflecting", 0.25, 250.0; T).x, 0.75; atol=1e-6)   # 2.75: off 1, off 0
            @test isapprox(step1("reflecting", 0.25, -300.0; T).x, 0.75; atol=1e-6)  # -2.75: off 0, off 1, off 0
            @test isapprox(step1("periodic", 0.25, 300.0; T).x, 0.25; atol=1e-6)
        end
        # the same for any step (no velocity), 2D and 3D; one reflection stays exact, as before
        @test DM.apply_boundary(mono(3.25), 1.0, "reflecting").x == 0.75
        e3 = DM.DiffusingEmitter3D{Float64}(0.5, -1.5, 2.25, 10.0, 0.0, 1, 1, 1, :monomer, nothing)
        @test coords(DM.apply_boundary(e3, 1.0, "reflecting")) == (0.5, 0.5, 0.25)
        @test DM.apply_boundary(mono(-0.1), 1.0, "reflecting").x == 0.1
        @test DM.apply_boundary(mono(1.1), 1.0, "reflecting").x == 2.0 - 1.1
        # a bound pair longer than the box: each end reflected on its own, folded the same way (-2.0 and 3.0)
        long = dm_config(box_size=1.0, boundary="reflecting", diff_dimer=0.0, diff_dimer_rot=0.0, k_off=0.0, d_dimer=5.0)
        lp = [DM.DiffusingEmitter2D{Float64}(0.2, 0.5, 10.0, 0.0, 1, 1, 1, :dimer, 2),
              DM.DiffusingEmitter2D{Float64}(0.8, 0.5, 10.0, 0.0, 1, 1, 2, :dimer, 1)]
        @test [e.x for e in DM.update_system(lp, long, 0.01)] == [0.0, 1.0]
        # a drift beyond Float32's range leaves a Float32 coordinate finite and inside the box (Codex review r2)
        for b in ("periodic", "reflecting")
            x = step1(b, 0.25, 1e41; T=Float32).x
            @test isfinite(x) && 0 <= x <= 1
        end
        e3f = DM.DiffusingEmitter3D{Float32}(0.25f0, 0.5f0, 0.5f0, 10.0, 0.0, 1, 1, 1, :monomer, nothing)
        p3 = dm_config(ndims=3, diff_monomer=0.0, box_size=1.0, r_react=0.0)
        @test all(c -> isfinite(c) && 0 <= c <= 1,
                  coords(only(DM.update_system([e3f], p3, 0.01; track_v=Dict(1 => (1e41, 0.0, -1e41))))))
        # the mean of two finite velocities never overflows, and a pair moving at it stays finite (Codex review r2)
        @test DM._pair_velocity((1e308, -1e308), (1e308, -1e308)) == (1e308, -1e308)
        @test DM._pair_velocity((1.0, 3.0), (2.0, -1.0)) == (1.5, 1.0)
        u = nextfloat(0.0)  # nor does it underflow (Codex review r3)
        @test DM._pair_velocity((u, 0.0), (u, 0.0)) == (u, 0.0)
        # the boundary is computed as on main for any coordinate type, BigFloat included (Codex review r3)
        eb = DM.DiffusingEmitter2D{BigFloat}(BigFloat(1) - BigFloat(2)^-60, BigFloat(0.5), 10.0, 0.0, 1, 1, 1, :monomer, nothing)
        @test DM.apply_boundary(eb, 1.0, "periodic").x == 1
        pr = [DM.DiffusingEmitter2D{Float64}(0.5, 0.5, 10.0, 0.0, 1, 1, 1, :dimer, 2),
              DM.DiffusingEmitter2D{Float64}(0.55, 0.5, 10.0, 0.0, 1, 1, 2, :dimer, 1)]
        pp = dm_config(box_size=1.0, diff_dimer=0.0, diff_dimer_rot=0.0, k_off=0.0, d_dimer=0.05)
        out = DM.update_system(pr, pp, 1e-308; track_v=Dict(1 => (1e308, 0.0), 2 => (1e308, 0.0)))
        @test length(out) == 2 && all(e -> all(c -> isfinite(c) && 0 <= c <= 1, coords(e)), out)
    end

    @testset "continuation" begin
        Random.seed!(9170)
        mix = [(0.5, 0.01), (0.5, 0.05)]
        p = dm_config(monomer_mobility=mix, monomer_velocity=[(0.5, 0.0), (0.0, -0.5)], t_max=0.1)
        smld, _ = simulate(p; override_count=100, γ=1000.0, track_velocity=Dict(7 => (0.0, 2.0)))
        ext = extract_end_state(smld)
        @test ext.metadata["monomer_velocity"] == smld.metadata["monomer_velocity"]
        # a saved velocity is a property of the molecule: kept under a config without velocities
        q = dm_config(monomer_mobility=mix, t_max=0.3)
        smld2, _ = @test_logs (:warn, r"keep their saved velocity") match_mode=:any simulate(q; starting_conditions=smld, γ=1000.0)
        @test smld2.metadata["monomer_velocity"] == smld.metadata["monomer_velocity"]
        steps = monomer_steps(smld2, q.box_size)
        @test isapprox(mean(d[2] for d in steps[7]), 0.02; atol=5 * sqrt(2 * 0.05 * 0.01 / length(steps[7])))
        # a run that set no velocity: its tracks keep their D and class and take this config's velocity for it
        src, _ = simulate(q; override_count=10, γ=1000.0)
        cont, _ = @test_logs simulate(p; starting_conditions=src)
        @test all(cont.metadata["monomer_velocity"][id] == p.monomer_velocity[k] for (id, k) in src.metadata["monomer_class"])
        # a changed velocity: kept, with one warning
        p2 = dm_config(monomer_mobility=mix, monomer_velocity=[(2.0, 0.0), (0.0, 2.0)], t_max=0.1)
        cont2, _ = @test_logs (:warn, r"keep their saved velocity") match_mode=:any simulate(p2; starting_conditions=smld)
        @test cont2.metadata["monomer_velocity"] == smld.metadata["monomer_velocity"]
        # the same config: kept, no warning; a molecule whose D was set with track_D keeps moving (review r1, M1)
        u = dm_config(monomer_velocity=[(0.5, 0.0)], t_max=0.3)
        s1, _ = simulate(u; override_count=5, γ=1000.0, track_D=Dict(1 => 0.0))
        s2, _ = @test_logs simulate(u; starting_conditions=s1)
        @test all(d -> all(isapprox.(d, (0.005, 0.0); atol=1e-12)), monomer_steps(s2, u.box_size)[1])
        # a tie (unknown provenance) carries no velocity
        tie = BasicSMLD(smld.emitters, smld.camera, smld.n_frames, smld.n_datasets,
                        merge(copy(smld.metadata), Dict{String,Any}("concatenated_from" => 2)))
        ext2 = @test_logs (:warn, r"concatenation") match_mode=:any extract_end_state(tie)
        @test !haskey(ext2.metadata, "monomer_velocity")
        # the identical config after a run that named a track: kept, no warning (review r2, S1)
        k = dm_config(t_max=0.3)
        k1, _ = simulate(k; override_count=5, γ=1000.0, track_velocity=Dict(1 => (0.5, 0.0)), track_D=Dict(1 => 0.0))
        @test k1.metadata["monomer_velocity_config"] == []
        k2, _ = @test_logs simulate(k; starting_conditions=k1)
        @test all(d -> all(isapprox.(d, (0.005, 0.0); atol=1e-12)), monomer_steps(k2, k.box_size)[1])
        # a saved velocity that does not fit is an ArgumentError naming its source
        md = copy(k1.metadata)
        md["monomer_velocity"] = merge(copy(k1.metadata["monomer_velocity"]), Dict(2 => (1.0, 0.0, 0.0)))
        bad = BasicSMLD(k1.emitters, k1.camera, k1.n_frames, k1.n_datasets, md)
        @test_throws r"saved velocity of track 2" simulate(k; starting_conditions=bad)
        # a saved D outside the model's range is an ArgumentError naming its source
        mdD = copy(k1.metadata)
        mdD["monomer_D"] = merge(copy(k1.metadata["monomer_D"]), Dict(2 => 1e60))
        @test_throws r"saved D of track 2" simulate(k; starting_conditions=BasicSMLD(k1.emitters, k1.camera,
                                                                                    k1.n_frames, k1.n_datasets, mdD))
        # a saved velocity that track_velocity replaces is not read (Codex review r3)
        md2 = copy(k1.metadata)
        md2["monomer_velocity"] = merge(copy(k1.metadata["monomer_velocity"]), Dict(2 => (1e308, 0.0)))
        huge = BasicSMLD(k1.emitters, k1.camera, k1.n_frames, k1.n_datasets, md2)
        k10 = dm_config(dt=10.0, t_max=20.0, camera_framerate=0.1, camera_exposure=10.0, k_off=0.0)
        @test_throws r"saved velocity of track 2" simulate(k10; starting_conditions=huge)
        s10, _ = simulate(k10; starting_conditions=huge, track_velocity=Dict(2 => (0.0, 0.0)))
        @test s10.metadata["monomer_velocity"][2] == (0.0, 0.0)
        # a kept class that is not one of this config's populations takes no velocity
        src3, _ = simulate(dm_config(monomer_mobility=[(0.3, 0.01), (0.3, 0.02), (0.4, 0.03)], t_max=0.02);
                           override_count=40, γ=1000.0)
        two = dm_config(monomer_mobility=[(0.5, 0.01), (0.5, 0.02)], monomer_velocity=[(1.0, 0.0), (0.0, 1.0)], t_max=0.02)
        c3, _ = @test_logs (:warn, r"monomer_mobility differs") match_mode=:any simulate(two; starting_conditions=src3)
        cls3, mv3 = src3.metadata["monomer_class"], c3.metadata["monomer_velocity"]
        @test any(==(3), values(cls3))
        @test all(mv3[id] == (cls3[id] == 3 ? (0.0, 0.0) : two.monomer_velocity[cls3[id]]) for id in keys(cls3))
        # a kept D that is no population's gives no class, so no population's velocity (Codex review r2: class 1
        # was invented, and with it population 1's velocity)
        s0, _ = simulate(dm_config(t_max=0.02); override_count=1, γ=1000.0, track_D=Dict(1 => 0.0))
        pops = dm_config(monomer_mobility=[(0.5, 0.1), (0.5, 0.2)], monomer_velocity=[(1.0, 0.0), (0.0, 1.0)], t_max=0.05)
        c0, _ = @test_logs (:warn, r"monomer_mobility differs") match_mode=:any simulate(pops; starting_conditions=s0)
        @test !haskey(c0.metadata["monomer_class"], 1) && c0.metadata["monomer_velocity"][1] == (0.0, 0.0)
        @test allequal(coords(e) for e in c0.emitters)  # D 0 and no velocity: it stays put
        c1, _ = @test_logs simulate(pops; starting_conditions=c0)
        @test !haskey(c1.metadata["monomer_class"], 1) && allequal(coords(e) for e in c1.emitters)
        # the two warnings count only tracks that keep their saved velocity or D: none when the keywords name them all
        # (Codex review r2)
        allv = Dict(id => (3.0, 0.0) for id in keys(smld.metadata["monomer_velocity"]))
        c4, _ = @test_logs simulate(p2; starting_conditions=smld, track_velocity=allv)
        @test all(==((3.0, 0.0)), values(c4.metadata["monomer_velocity"]))
        q5 = dm_config(monomer_mobility=[(1.0, 0.2)], t_max=0.1)
        allD = Dict(id => 0.1 for id in keys(smld.metadata["monomer_D"]))
        @test_logs simulate(q5; starting_conditions=smld, track_velocity=allv, track_D=allD)
        @test_logs (:warn, r"keep their saved D") simulate(q5; starting_conditions=smld, track_velocity=allv)
        @test_logs (:warn, r"keep their saved velocity") simulate(q5; starting_conditions=smld, track_D=allD)
    end

    @testset "binding and unbinding of directed molecules" begin
        mono(id, x) = DM.DiffusingEmitter2D{Float64}(x, 1.0, 10.0, 0.0, 1, 1, id, :monomer, nothing)
        bound(id, x, partner) = DM.DiffusingEmitter2D{Float64}(x, 1.0, 10.0, 0.0, 1, 1, id, :dimer, partner)
        cfg(pm; kw...) = dm_config(; r_react=0.01, d_dimer=0.02, k_off=0.0, t_max=0.03, pair_mobility=pm, kw...)
        at(s, id, t) = only(e for e in s.emitters if e.track_id == id && isapprox(e.timestamp, t; atol=1e-12))
        sim(p, start; kw...) = (Random.seed!(9180); simulate(p; starting_conditions=start, γ=1000.0, kw...)[1])
        D0 = Dict(1 => 0.0, 2 => 0.0)
        # two directed molecules with D 0 form a pair: neither is immobile, so under :min they are placed about
        # their midpoint (not anchored), with no drift on the formation step; under :fixed there is no warning
        for pm in (:min, :fixed)
            s = @test_logs sim(cfg(pm), [mono(1, 1.0), mono(2, 1.004)]; track_D=D0,
                               track_velocity=Dict(1 => (1.0, 0.0), 2 => (0.0, 0.5)))
            a, b = at(s, 1, 0.01), at(s, 2, 0.01)
            @test a.state == :dimer && isapprox(a.x, 0.992; atol=1e-12) && isapprox(b.x, 1.012; atol=1e-12)
            @test a.y == 1.0 && b.y == 1.0
        end
        # a directed D 0 molecule meets an immobile one under :min: the immobile one anchors
        s = sim(cfg(:min), [mono(1, 1.0), mono(2, 1.004)]; track_D=D0, track_velocity=Dict(1 => (1.0, 0.0)))
        @test at(s, 2, 0.01).x == 1.004 && isapprox(at(s, 1, 0.01).x, 0.984; atol=1e-12)
        @test at(s, 2, 0.02).x == 1.004
        # a split closer than r_react under :min: the immobile partner stays, the directed one is placed and drifts
        sp(pm) = cfg(pm; d_dimer=0.004, k_off=100.0)
        s = sim(sp(:min), [bound(1, 1.0, 2), bound(2, 1.004, 1)]; track_D=D0, track_velocity=Dict(1 => (-1.0, 0.0)))
        @test at(s, 1, 0.01).state == :monomer && at(s, 2, 0.01).x == 1.004
        @test at(s, 1, 0.01).x < 1.004 - 0.01 - 0.01 + 1e-9
        # a split under :fixed of two directed D 0 molecules moves both, with no warning
        s = @test_logs sim(sp(:fixed), [bound(1, 1.0, 2), bound(2, 1.004, 1)]; track_D=D0,
                           track_velocity=Dict(1 => (-1.0, 0.0), 2 => (1.0, 0.0)))
        @test at(s, 1, 0.01).state == :monomer
        # the dissociation step drifts each partner (a pair wider than r_react is not re-placed)
        s = sim(cfg(:fixed; d_dimer=0.05, k_off=100.0), [bound(1, 1.0, 2), bound(2, 1.05, 1)]; track_D=D0,
                track_velocity=Dict(1 => (1.0, 0.0), 2 => (0.0, 1.0)))
        @test isapprox(at(s, 1, 0.01).x, 1.01; atol=1e-12) && at(s, 1, 0.01).y == 1.0
        @test at(s, 2, 0.01).x == 1.05 && isapprox(at(s, 2, 0.01).y, 1.01; atol=1e-12)
    end

    @testset "partners of different coordinate types each keep their own (Codex review r7)" begin
        # main threw MethodError when a Float64 and a Float32 molecule formed or kept a pair, and dissociation
        # converted the partner to the first molecule's type
        e(T, id, x; state=:monomer, partner=nothing, ph=10.0) =
            DM.DiffusingEmitter2D{T}(T(x), T(0.5), T(ph), T(0), 1, 1, id, state, partner)
        types(out) = Dict(o.track_id => typeof(o.x) for o in out)
        xs(out) = Dict(o.track_id => Float64(o.x) for o in out)
        D0 = Dict(1 => 0.0, 2 => 0.0)  # both immobile (pair_mobility = :min needs diff_monomer > 0)
        for b in ("periodic", "reflecting"), pm in (:fixed, :min), (T1, T2) in ((Float64, Float32), (Float32, Float64))
            # formation (both immobile: :min anchors track 1), mid-box and at a wall, then a bound step
            for x1 in (0.5, 0.0)
                p = dm_config(box_size=1.0, boundary=b, r_react=0.01, d_dimer=0.02, k_off=0.0, pair_mobility=pm)
                out = DM.update_system(DM.AbstractDiffusingEmitter[e(T1, 1, x1), e(T2, 2, x1 + 0.005)], p, 0.01;
                                       track_D=D0)
                @test types(out) == Dict(1 => T1, 2 => T2)
                @test all(o -> o.state == :dimer, out)
                d = xs(out)[2] - xs(out)[1]
                @test isapprox(b == "periodic" ? d - round(d) : d, 0.02; atol=1e-6)
                out = DM.update_system(out, p, 0.01; track_D=D0)
                @test types(out) == Dict(1 => T1, 2 => T2) && all(o -> o.state == :dimer, out)
            end
            # dissociation of a pair closer than r_react: each keeps its type, and they end r_react apart in
            # the stored coordinates (the margin is the coarser type's)
            p = dm_config(box_size=1.0, boundary=b, r_react=0.01, d_dimer=0.004, k_off=Inf, pair_mobility=pm)
            for k in 1:20
                x2 = 0.5 + 0.004 + k * 1.37e-5
                out = DM.update_system(DM.AbstractDiffusingEmitter[e(T1, 1, 0.5; state=:dimer, partner=2),
                                                                 e(T2, 2, x2; state=:dimer, partner=1)], p, 0.01;
                                       track_D=D0)
                @test types(out) == Dict(1 => T1, 2 => T2) && all(o -> o.state == :monomer, out)
                @test abs(xs(out)[2] - xs(out)[1]) >= 0.01
            end
        end
        # a bound step that takes an end past a reflecting wall: the pair is re-centered off the wall, each end in
        # its own type
        for (T1, T2) in ((Float64, Float32), (Float32, Float64))
            p = dm_config(box_size=1.0, boundary="reflecting", d_dimer=0.02, k_off=0.0)
            o1, o2 = DM._move_pair(e(T1, 1, 0.5; state=:dimer, partner=2), e(T2, 2, 0.5; state=:dimer, partner=1),
                                   (0.005, 0.5), (0.01, 0.0), 0.02, p)
            @test o1 isa DM.DiffusingEmitter2D{T1} && o2 isa DM.DiffusingEmitter2D{T2}
            # the center 0.005 folds off 0.01 (the wall moved in by the half-bond) to 0.015
            @test isapprox(Float64(o1.x), 0.005; atol=1e-7) && isapprox(Float64(o2.x), 0.025; atol=1e-7)
        end
        # Codex's case: a Float64 partner keeps its coordinates and photons when its Float32 partner dissociates
        p = dm_config(box_size=1.0, diff_monomer=0.0, r_react=0.0, k_off=Inf)
        out = DM.update_system(DM.AbstractDiffusingEmitter[e(Float32, 1, 0.5; state=:dimer, partner=2),
                                                         e(Float64, 2, 0.50000001; state=:dimer, partner=1, ph=1e40)], p, 0.01)
        o2 = only(o for o in out if o.track_id == 2)
        @test o2 isa DM.DiffusingEmitter2D{Float64} && o2.x == 0.50000001 && o2.photons == 1e40
        e3(T, id, x, partner) = DM.DiffusingEmitter3D{T}(T(x), T(0.5), T(0.5), T(10), T(0), 1, 1, id, :dimer, partner)
        out = DM.update_system(DM.AbstractDiffusingEmitter[e3(Float32, 1, 0.5, 2), e3(Float64, 2, 0.50000001, 1)],
                               dm_config(ndims=3, box_size=1.0, diff_monomer=0.0, r_react=0.0, k_off=Inf), 0.01)
        @test types(out) == Dict(1 => Float32, 2 => Float64)
        @test only(o for o in out if o.track_id == 2).x == 0.50000001
        # a whole run: pairs form and split, and every record keeps its molecule's type
        Random.seed!(9190)
        s, _ = simulate(DiffusionSMLMConfig(box_size=1.0, t_max=0.1, r_react=0.01);
                        starting_conditions=[e(Float64, 1, 0.5), e(Float32, 2, 0.501)], γ=1000.0)
        @test any(r -> r.state == :dimer, s.emitters)
        @test all(r -> typeof(r.x) == (r.track_id == 1 ? Float64 : Float32), s.emitters)
    end

    @testset "validation" begin
        @test_throws ArgumentError DiffusionSMLMConfig(monomer_velocity=[(1.0, 0.0), (0.0, 1.0)])
        @test_throws ArgumentError DiffusionSMLMConfig(monomer_velocity=[(1.0, 0.0, 0.0)])
        @test_throws ArgumentError DiffusionSMLMConfig(ndims=3, monomer_velocity=[(1.0, 0.0)])
        @test_throws ArgumentError DiffusionSMLMConfig(monomer_velocity=[(NaN, 0.0)])
        @test_throws ArgumentError DiffusionSMLMConfig(monomer_mobility=[(0.5, 0.1), (0.5, 0.2)],
                                                       monomer_velocity=[(1.0, 0.0)])
        p = dm_config(t_max=0.02)
        @test_throws ArgumentError simulate(p; override_count=3, track_velocity=Dict(4 => (1.0, 0.0)))
        @test_throws ArgumentError simulate(p; override_count=3, track_velocity=Dict(1 => (1.0, 0.0, 0.0)))
        @test_throws ArgumentError simulate(p; override_count=3, track_velocity=Dict(1 => (Inf, 0.0)))
        @test_throws ArgumentError simulate(p; override_count=3, track_D=Dict(1 => -0.1))
        @test_throws ArgumentError simulate(p; override_count=3, track_D=Dict(1 => NaN))
        @test_throws ArgumentError simulate(p; override_count=3, track_D=Dict(9 => 0.1))
        @test_throws ArgumentError simulate(p; override_count=3, track_velocity=Dict(1.5 => (1.0, 0.0)))
        @test_throws ArgumentError simulate(p; override_count=3, track_velocity=Dict(true => (1.0, 0.0)))
        @test_throws ArgumentError simulate(p; override_count=3, track_D=Dict(1e30 => 0.1))
        @test_throws ArgumentError simulate(p; override_count=3, track_velocity=Dict(1 => nothing))
        @test_throws ArgumentError simulate(p; override_count=3, track_D=Dict("1" => 0.1))
        # every scalar field: NaN throws, Inf throws except for k_off and r_react, and so does a value that converts
        # to Inf (Codex review r2)
        for f in (:density, :box_size, :diff_monomer, :diff_dimer, :diff_dimer_rot, :k_off, :r_react, :d_dimer, :dt,
                  :t_max, :camera_framerate, :camera_exposure)
            @test_throws ArgumentError DiffusionSMLMConfig(; f => NaN)
            f in (:k_off, :r_react) || @test_throws ArgumentError DiffusionSMLMConfig(; f => Inf)
        end
        @test DiffusionSMLMConfig(k_off=Inf, r_react=Inf).k_off == Inf
        @test_throws ArgumentError DiffusionSMLMConfig(diff_monomer=big"1e400")
        @test_throws ArgumentError DiffusionSMLMConfig(monomer_mobility=[(0.5, big"1e400"), (0.5, 0.1)])
        @test_throws ArgumentError simulate(p; override_count=3, track_D=Dict(1 => big"1e400"))
        # a malformed velocity entry is an ArgumentError, not a MethodError (Codex review r2)
        @test_throws ArgumentError DiffusionSMLMConfig(monomer_velocity=[nothing])
        @test_throws ArgumentError DiffusionSMLMConfig(monomer_velocity=[(1.0, "a")])
        @test_throws ArgumentError DiffusionSMLMConfig(monomer_velocity=(1.0, 0.0))  # one tuple, not a vector of them
        # the model's range: 0, or a magnitude from 1e-50 to 1e50, for every number but r_react (Codex review r4:
        # sqrt(2·D·dt) overflowed at D = 1e308 and underflowed at D = 1e-300 within the old wording)
        for f in (:density, :box_size, :diff_monomer, :diff_dimer, :diff_dimer_rot, :k_off, :d_dimer, :dt, :t_max,
                  :camera_framerate, :camera_exposure), x in (1e60, 1e-60)
            @test_throws ArgumentError DiffusionSMLMConfig(; f => x)
        end
        @test DiffusionSMLMConfig(r_react=1e-201).r_react == 1e-201  # compared with a scaled distance only
        @test DiffusionSMLMConfig(diff_monomer=0.0, diff_dimer=0.0, diff_dimer_rot=0.0, k_off=0.0).k_off == 0.0
        @test_throws ArgumentError DiffusionSMLMConfig(monomer_velocity=[(1e-60, 0.0)])
        @test_throws ArgumentError DiffusionSMLMConfig(monomer_velocity=[(1e308, 0.0)])
        @test_throws ArgumentError DiffusionSMLMConfig(monomer_mobility=[(0.5, 1e-60), (0.5, 0.1)])
        @test_throws ArgumentError DM._check_track_velocity(Dict(1 => (1e60, 0.0)), 2)
        @test_throws ArgumentError simulate(p; override_count=3, track_D=Dict(1 => 1e60))
        # (Codex review r5) checked as given, before conversion; mixture fractions too
        @test_throws ArgumentError DiffusionSMLMConfig(monomer_mobility=[(1e-60, 0.0), (1.0, 0.1)])
        @test_throws ArgumentError DiffusionSMLMConfig(monomer_velocity=[(big"1e-400", 0)])
        @test_throws ArgumentError simulate(p; override_count=3, track_D=Dict(1 => -big"1e-400"))
        # simulate checks the whole config again, on its fields as they are
        for change! in (q -> q.diff_monomer = 1e308, q -> q.monomer_mobility = [(1.0, 1e308)], q -> q.dt = 1e308)
            q = dm_config(t_max=0.02)
            change!(q)
            @test_throws ArgumentError simulate(q; override_count=1)
        end
        # γ in the model's range, so γ·dt stays normal
        @test_throws ArgumentError simulate(p; override_count=1, γ=1e308)
        @test_throws ArgumentError simulate(p; override_count=1, γ=1e-300)
        # an exposure far longer than the frame is capped before it is counted; a run beyond 2^53 steps throws
        long_exposure = dm_config(dt=0.1, camera_framerate=10.0, camera_exposure=1e50, t_max=0.1)
        s_long, _ = @test_logs (:warn, r"capping the exposure") match_mode=:any simulate(long_exposure; override_count=1)
        @test length(s_long.emitters) == 1
        @test_throws ArgumentError simulate(dm_config(dt=1e-50, t_max=1.0, camera_framerate=1.0, camera_exposure=1.0);
                                           override_count=1)
        @test_throws ArgumentError simulate(dm_config(dt=1.0, t_max=1e50, camera_framerate=1.0, camera_exposure=1.0);
                                           override_count=1)
        # Float32 starting conditions: the box, a record's photons and the run's duration must fit Float32
        f32 = [DM.DiffusingEmitter2D{Float32}(0.5f0, 0.5f0, 10f0, 0f0, 1, 1, 1, :monomer, nothing)]
        long32 = dm_config(box_size=1.0, diff_monomer=0.0, r_react=0.0, k_off=0.0, dt=1e40, camera_exposure=1e40,
                           t_max=2e40, camera_framerate=1e-40)
        @test_throws ArgumentError simulate(long32; starting_conditions=f32)
        # (Codex review r6) every starting molecule's type is checked, in either order
        e64 = DM.DiffusingEmitter2D{Float64}(0.5, 0.5, 10.0, 0.0, 1, 1, 1, :monomer, nothing)
        e32 = DM.DiffusingEmitter2D{Float32}(0.5f0, 0.5f0, 10f0, 0f0, 1, 1, 2, :monomer, nothing)
        for mixed in (DM.AbstractDiffusingEmitter[e64, e32], DM.AbstractDiffusingEmitter[e32, e64])
            @test_throws ArgumentError simulate(long32; starting_conditions=mixed, camera=IdealCamera(1:1, 1:1, 1.0))
        end
        # the deprecated photons keyword stays unchecked, Inf included, as in 0.7
        s_inf, _ = @test_logs (:warn, r"photons keyword is deprecated") match_mode=:any simulate(
            DiffusionSMLMConfig(t_max=0.1); override_count=1, photons=Inf)
        @test all(e -> e.photons == Inf, s_inf.emitters)
        # a velocity that is the smallest integer of its type is in the range (abs would overflow)
        @test DiffusionSMLMConfig(monomer_velocity=[(Int8(-128), 0)]).monomer_velocity == [(-128.0, 0.0)]
        # counts compared exactly: 3 frames × 3002399751580331 steps is 2^53 + 1, which Float64 rounds to 2^53 (checked
        # on the frame count directly, so a regression fails here instead of starting that run)
        dtx = 1 / 3002399751580331
        cx = dm_config(dt=dtx, camera_exposure=dtx, camera_framerate=1.0, t_max=3.0)
        @test DM.substeps_per_frame(cx) == (1, 3002399751580331)
        @test_throws ArgumentError DM._frame_count(cx, 3002399751580331)
        @test DM._frame_count(cx, 3002399751580330) == 3  # 2^53 - 2 steps
        # a default camera or a molecule count beyond 2^53 is an ArgumentError, not an InexactError
        @test_throws ArgumentError simulate(DiffusionSMLMConfig(box_size=1e20, t_max=0.1); override_count=1)
        @test_throws ArgumentError simulate(DiffusionSMLMConfig(density=1e20, box_size=1.0, t_max=0.1))
        # the config is mutable: simulate checks its velocities again
        p.monomer_velocity = [(1.0, 0.0), (0.0, 1.0)]
        @test_throws ArgumentError simulate(p; override_count=3)
        # positional construction: the 16-argument form (before monomer_velocity) and the full one
        args = (1.0, 10.0, 0.1, 0.05, 0.5, 0.2, 0.01, 0.05, 0.01, 10.0, 2, "periodic", 10.0, 0.1,
                Tuple{Float64,Float64}[], :fixed)
        @test isempty(DiffusionSMLMConfig(args...).monomer_velocity)
        @test DiffusionSMLMConfig(args..., [(0.1, 0.2)]).monomer_velocity == [(0.1, 0.2)]
    end
end
