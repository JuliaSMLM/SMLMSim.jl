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
            # the config's one velocity is not saved per track (as diff_monomer); the key says the run has one
            @test haskey(smld.metadata, "monomer_velocity") && isempty(smld.metadata["monomer_velocity"])
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
        @test sort!(collect(keys(smld.metadata["monomer_velocity"]))) == 1:200
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

    @testset "reflecting wall" begin
        Random.seed!(9160)
        p = dm_config(diff_monomer=0.0, box_size=2.0, boundary="reflecting", monomer_velocity=[(2.0, 0.0)], t_max=1.5)
        smld, _ = simulate(p; override_count=20, γ=1000.0)
        @test all(e -> 0 <= e.x <= p.box_size, smld.emitters)
        last_x = [e.x for e in smld.emitters if e.frame == smld.n_frames]
        @test all(x -> x >= p.box_size - 2.0 * p.dt - 1e-12, last_x)  # held within one step of the wall
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
        smld2, _ = simulate(q; starting_conditions=smld, γ=1000.0)
        @test smld2.metadata["monomer_velocity"] == smld.metadata["monomer_velocity"]
        steps = monomer_steps(smld2, q.box_size)
        @test isapprox(mean(d[2] for d in steps[7]), 0.02; atol=5 * sqrt(2 * 0.05 * 0.01 / length(steps[7])))
        # a run that set no velocity: its tracks keep their D and stay without drift, with a warning
        src, _ = simulate(q; override_count=10, γ=1000.0)
        cont = @test_logs (:warn, r"keep the source run's velocity") match_mode=:any simulate(p; starting_conditions=src)
        @test all(v -> v == (0.0, 0.0), values(cont[1].metadata["monomer_velocity"]))
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
